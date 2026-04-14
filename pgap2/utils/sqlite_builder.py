"""Build a self-contained SQLite database from PGAP2 output files.

This is the **data-producer** side: PGAP2 runs its analysis pipeline and
then serialises the results into a portable SQLite ``.db`` file.  The
companion web application (PGAP2-Web) is a **data-consumer** that reads
these SQLite files and registers lightweight metadata in MySQL.

Usage (programmatic)::

    from pgap2.utils.sqlite_builder import build_sqlite
    info = build_sqlite(data_dir="/path/to/pgap2/output",
                        output_path="/path/to/pgap2.db",
                        project_name="My Project",
                        species="E. coli")

Usage (CLI)::

    pgap2 sqlite -i /path/to/output -o /path/to/output --name "My Project"
"""

from __future__ import annotations

import csv
import json
import re
import zlib
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from loguru import logger
from sqlalchemy import create_engine, event, select
from sqlalchemy.orm import Session, sessionmaker

from pgap2.models import (
    Base, Project, Strain, GeneCluster, PanGroup, PAV,
    GeneAnnotation,
    Rarefaction, PanGroupStat, CurveFit, ClusterStrainFreq,
    ParalogStat, StatAttrBin,
    PhylogenyTree, BAPSCluster, TajimasD,
    PrepStat, GeneCodeUsage,
    ProjectFile,
)

# Map the group labels used in detail.tsv to PanGroup enum values
GROUP_MAP = {
    "Strict_core": PanGroup.STRICT_CORE,
    "Core": PanGroup.CORE,
    "Soft_core": PanGroup.SOFT_CORE,
    "Shell": PanGroup.SHELL,
    "Cloud": PanGroup.CLOUD,
}


# ═══════════════════════════════════════════════════════
# Public API
# ═══════════════════════════════════════════════════════

def build_sqlite(
    data_dir: str,
    output_path: str | None = None,
    project_name: str = "PGAP2 Project",
    species: str = "Unknown",
    description: str = "",
) -> dict:
    """Parse PGAP2 output files and write everything into a SQLite DB.

    Parameters
    ----------
    data_dir : str
        Path to the PGAP2 output directory (contains *.pav, *.detail.tsv, …).
    output_path : str, optional
        Path for the resulting SQLite file.  Defaults to
        ``<data_dir>/pgap2.db``.
    project_name : str
        Human-readable project name.
    species : str
        Species name.
    description : str
        Optional description text.

    Returns
    -------
    dict
        Keys: sqlite_path, project_id, num_strains, num_clusters,
              data_dir, strain_order, strain_info.
    """
    data_path = Path(data_dir).resolve()
    if not data_path.is_dir():
        raise FileNotFoundError(f"Data directory not found: {data_dir}")

    if output_path is None:
        output_path = str(data_path / "pgap2.db")
    sqlite_file = Path(output_path)
    sqlite_file.parent.mkdir(parents=True, exist_ok=True)

    # Remove stale file if present
    if sqlite_file.exists():
        sqlite_file.unlink()

    logger.info("Creating database: {}", output_path)

    engine = create_engine(f"sqlite:///{output_path}", echo=False)

    @event.listens_for(engine, "connect")
    def _set_sqlite_pragmas(dbapi_conn, _record):
        cursor = dbapi_conn.cursor()
        cursor.execute("PRAGMA journal_mode=WAL")
        cursor.execute("PRAGMA foreign_keys=OFF")  # OFF during bulk insert
        cursor.close()

    Base.metadata.create_all(engine)

    SessionLocal = sessionmaker(bind=engine, expire_on_commit=False)

    with SessionLocal() as db:
        info = _import_all(db, data_path, project_name, species, description)
        db.commit()

    engine.dispose()
    info["sqlite_path"] = str(sqlite_file)

    db_size = sqlite_file.stat().st_size
    size_mb = db_size / (1024 * 1024)
    logger.success(
        "Done! {} strains, {} clusters -> {:.1f} MB",
        info["num_strains"], info["num_clusters"], size_mb,
    )
    return info


# ═══════════════════════════════════════════════════════
# Internal: orchestrate all import steps
# ═══════════════════════════════════════════════════════

def _import_all(
    db: Session,
    data_path: Path,
    project_name: str,
    species: str,
    description: str,
) -> dict:
    """Run every import phase inside a single session."""

    # 1. Create project
    logger.info("[1/9] Creating project ...")
    project = Project(
        name=project_name,
        species=species,
        description=description,
        data_dir=str(data_path),
    )
    db.add(project)
    db.flush()
    logger.debug("Created project {} (id={})", project_name, project.id)

    # 2. Parse detail TSV → strains + clusters + PAV
    logger.info("[2/9] Importing strains, clusters & PAV ...")
    detail_file = data_path / "pgap2.partition.gene_content.detail.tsv"
    pav_file = data_path / "pgap2.partition.gene_content.pav"

    strain_map: Dict[str, Strain] = {}
    cluster_map: Dict[str, GeneCluster] = {}
    strain_names: List[str] = []

    if detail_file.exists():
        strain_map, cluster_map, strain_names = _import_detail_tsv(
            db, project.id, detail_file,
            pav_file if pav_file.exists() else None,
        )
    elif pav_file.exists():
        strain_map, cluster_map, strain_names = _import_pav_only(
            db, project.id, pav_file)

    project.num_strains = len(strain_map)
    project.num_clusters = len(cluster_map)
    project.strain_order = ",".join(strain_names) if strain_names else None
    logger.info("       {} strains, {} clusters loaded",
                len(strain_map), len(cluster_map))

    # 3. Summary statistics fallback
    logger.info("[3/9] Importing statistics ...")
    pan_group_file = _find_file(data_path,
                                "postprocess/postprocess.pan_group_stat.tsv",
                                "postprocess.pan_group_stat.tsv")
    summary_file = data_path / "pgap2.partition.summary_statistics.txt"
    if not pan_group_file and summary_file.exists():
        _import_summary_statistics(db, project.id, summary_file)

    # 4. Gene annotations
    logger.info("[4/9] Importing gene annotations ...")
    annot_file = data_path / "total.involved_annot.tsv"
    if annot_file.exists():
        _import_annotations(db, project.id, annot_file,
                            strain_map, cluster_map)

    # 5. Post-processing results
    logger.info("[5/9] Importing post-processing results ...")
    _import_post_results(db, project.id, data_path)

    # 5b. Cluster-strain frequency histogram
    freq_file = _find_file(data_path,
                           "postprocess/postprocess.clust_strain_freq.tsv",
                           "postprocess.clust_strain_freq.tsv")
    if not freq_file:
        _compute_freq_histogram(db, project.id, cluster_map)

    # 6. Trees
    logger.info("[6/9] Importing phylogeny trees ...")
    _import_trees(db, project.id, data_path)

    # 7. BAPS
    logger.info("[7/9] Importing BAPS clusters ...")
    _import_baps(db, project.id, data_path, strain_map)

    # 8. Preprocess stats
    logger.info("[8/9] Importing preprocess statistics ...")
    _import_prep_stats(db, project.id, data_path, strain_map)
    _import_gene_code(db, project.id, data_path)

    # 9. Embed disk files (MSA, GML, layout JSON …) as BLOBs
    logger.info("[9/9] Embedding files (MSA, graph, layout) ...")
    _embed_files(db, project.id, data_path)

    db.flush()

    # Collect strain info for downstream use
    strain_info = [
        (s.name, s.assembly_accession) for s in strain_map.values()
    ]

    logger.debug("Import complete: {} strains, {} clusters",
                 project.num_strains, project.num_clusters)

    return {
        "project_id": project.id,
        "num_strains": project.num_strains,
        "num_clusters": project.num_clusters,
        "data_dir": str(data_path),
        "strain_order": project.strain_order,
        "strain_info": strain_info,
    }


# ═══════════════════════════════════════════════════════
# Detail TSV + PAV parsing
# ═══════════════════════════════════════════════════════

def _import_detail_tsv(
    db: Session,
    project_id: int,
    detail_file: Path,
    pav_file: Optional[Path],
) -> Tuple[Dict[str, Strain], Dict[str, GeneCluster], List[str]]:
    """Parse gene_content.detail.tsv → strains, clusters, pav_vector."""
    strain_map: Dict[str, Strain] = {}
    cluster_map: Dict[str, GeneCluster] = {}

    pav_matrix: Dict[str, Dict[str, int]] = {}
    if pav_file:
        pav_matrix = _read_pav_matrix(pav_file)

    with open(detail_file, "r") as f:
        header = f.readline().strip()
        parts = header.split("\t")

        # Discover strain names (from PAV header or CSV header)
        strain_names: List[str] = []
        if pav_file:
            with open(pav_file, "r") as pf:
                pav_header = pf.readline().strip().split("\t")
                strain_names = pav_header[1:]
        else:
            csv_file = detail_file.parent / "pgap2.partition.gene_content.csv"
            if csv_file.exists():
                with open(csv_file, "r") as cf:
                    csv_header = cf.readline().strip().split(",")
                    strain_names = csv_header[1:]

        # Create strain records
        for sname in strain_names:
            strain = Strain(project_id=project_id, name=sname)
            if sname.startswith("GCF_") or sname.startswith("GCA_"):
                strain.assembly_accession = sname
            db.add(strain)
            strain_map[sname] = strain
        db.flush()

        # Parse cluster rows
        cluster_batch: List[GeneCluster] = []
        batch_size = 1000

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 13:
                continue

            cluster_id = parts[0]
            product = parts[2] if parts[2] else None
            group_str = parts[3]
            involved_strain = int(parts[9]) if parts[9] else 0
            involved_gene = int(parts[11]) if parts[11] else 0

            group = GROUP_MAP.get(group_str, PanGroup.CLOUD)
            frequency = (involved_strain / len(strain_names)
                         if strain_names else 0)

            avg_length = None
            try:
                avg_length = float(parts[6]) if parts[6] else None
            except ValueError:
                pass

            pav_vector = None
            if pav_matrix and cluster_id in pav_matrix:
                strain_vals = pav_matrix[cluster_id]
                pav_vector = ",".join(
                    str(strain_vals.get(sname, 0)) for sname in strain_names)

            cluster = GeneCluster(
                project_id=project_id,
                cluster_id=cluster_id,
                annotation=product,
                num_strains=involved_strain,
                num_genes=involved_gene,
                avg_length=avg_length,
                group=group,
                frequency=frequency,
                pav_vector=pav_vector,
            )
            cluster_batch.append(cluster)
            cluster_map[cluster_id] = cluster

            if len(cluster_batch) >= batch_size:
                db.add_all(cluster_batch)
                db.flush()
                cluster_batch = []

        if cluster_batch:
            db.add_all(cluster_batch)
            db.flush()

    # Create PAV rows
    if pav_matrix and strain_map and cluster_map:
        _create_pav_rows(db, project_id, pav_matrix, strain_map, cluster_map)

    return strain_map, cluster_map, strain_names


def _read_pav_matrix(pav_file: Path) -> Dict[str, Dict[str, int]]:
    """Read PAV binary matrix → {cluster_id: {strain: copy_number}}."""
    result: Dict[str, Dict[str, int]] = {}
    with open(pav_file, "r") as f:
        header = f.readline().strip().split("\t")
        strain_names = header[1:]
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            cluster_id = parts[0]
            vals = {}
            for i, sname in enumerate(strain_names):
                try:
                    vals[sname] = int(parts[i + 1])
                except (IndexError, ValueError):
                    vals[sname] = 0
            result[cluster_id] = vals
    return result


def _import_pav_only(
    db: Session, project_id: int, pav_file: Path,
) -> Tuple[Dict[str, Strain], Dict[str, GeneCluster], List[str]]:
    """Fallback: build strains/clusters from PAV matrix alone."""
    strain_map: Dict[str, Strain] = {}
    cluster_map: Dict[str, GeneCluster] = {}

    pav_matrix = _read_pav_matrix(pav_file)
    if not pav_matrix:
        return strain_map, cluster_map, []

    strain_names = list(next(iter(pav_matrix.values())).keys())
    for sname in strain_names:
        strain = Strain(project_id=project_id, name=sname)
        if sname.startswith("GCF_") or sname.startswith("GCA_"):
            strain.assembly_accession = sname
        db.add(strain)
        strain_map[sname] = strain
    db.flush()

    cluster_batch: List[GeneCluster] = []
    for cluster_id, strain_vals in pav_matrix.items():
        num_present = sum(1 for v in strain_vals.values() if v > 0)
        freq = num_present / len(strain_names) if strain_names else 0

        if freq >= 0.99:
            group = PanGroup.CORE
        elif freq >= 0.95:
            group = PanGroup.SOFT_CORE
        elif freq >= 0.15:
            group = PanGroup.SHELL
        else:
            group = PanGroup.CLOUD

        pav_vector = ",".join(
            str(strain_vals.get(sname, 0)) for sname in strain_names)

        cluster = GeneCluster(
            project_id=project_id,
            cluster_id=cluster_id,
            num_strains=num_present,
            num_genes=sum(strain_vals.values()),
            group=group,
            frequency=freq,
            pav_vector=pav_vector,
        )
        cluster_batch.append(cluster)
        cluster_map[cluster_id] = cluster

        if len(cluster_batch) >= 1000:
            db.add_all(cluster_batch)
            db.flush()
            cluster_batch = []

    if cluster_batch:
        db.add_all(cluster_batch)
        db.flush()

    _create_pav_rows(db, project_id, pav_matrix, strain_map, cluster_map)
    return strain_map, cluster_map, strain_names


def _create_pav_rows(
    db: Session,
    project_id: int,
    pav_matrix: Dict[str, Dict[str, int]],
    strain_map: Dict[str, Strain],
    cluster_map: Dict[str, GeneCluster],
):
    """Bulk-insert PAV rows from pav_matrix."""
    pav_batch: List[PAV] = []
    batch_size = 5000
    for cid, strain_vals in pav_matrix.items():
        cluster = cluster_map.get(cid)
        if not cluster:
            continue
        for sname, copies in strain_vals.items():
            if copies <= 0:
                continue
            strain = strain_map.get(sname)
            if not strain:
                continue
            pav_batch.append(PAV(
                project_id=project_id,
                gene_cluster_id=cluster.id,
                strain_id=strain.id,
                copy_number=copies,
            ))
            if len(pav_batch) >= batch_size:
                db.add_all(pav_batch)
                db.flush()
                pav_batch = []
    if pav_batch:
        db.add_all(pav_batch)
        db.flush()
    logger.debug("Created PAV rows for {} clusters × {} strains",
                 len(cluster_map), len(strain_map))


# ═══════════════════════════════════════════════════════
# Summary statistics
# ═══════════════════════════════════════════════════════

def _import_summary_statistics(
    db: Session, project_id: int, summary_file: Path,
):
    """Parse summary_statistics.txt → PanGroupStat rows."""
    group_name_map = {
        "Strict core genes": "strict_core",
        "Core genes": "core",
        "Soft core genes": "soft_core",
        "Shell genes": "shell",
        "Cloud genes": "cloud",
    }
    group_stats: Dict[str, Tuple[int, float]] = {}

    with open(summary_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            category = parts[0].strip()
            try:
                count = int(parts[2].strip())
            except (ValueError, IndexError):
                continue
            pct_str = parts[3].strip().rstrip("%")
            pct = float(pct_str) if pct_str else 0

            group = group_name_map.get(category)
            if group and group != "total":
                if group in group_stats:
                    old_count, old_pct = group_stats[group]
                    group_stats[group] = (old_count + count, old_pct + pct)
                else:
                    group_stats[group] = (count, pct)

    for group, (count, pct) in group_stats.items():
        db.add(PanGroupStat(
            project_id=project_id,
            group=group,
            num_clusters=count,
            num_genes=count,
            percentage=pct,
        ))
    db.flush()


# ═══════════════════════════════════════════════════════
# Gene annotations
# ═══════════════════════════════════════════════════════

def _import_annotations(
    db: Session,
    project_id: int,
    annot_file: Path,
    strain_map: Dict[str, Strain],
    cluster_map: Dict[str, GeneCluster],
):
    """Stream-import gene annotations from total.involved_annot.tsv."""
    # Build gene_id → cluster_id index from the CSV
    csv_file = annot_file.parent / "pgap2.partition.gene_content.csv"
    gene_to_cluster: Dict[str, str] = {}
    if csv_file.exists():
        with open(csv_file, "r") as f:
            f.readline()  # header
            for line in f:
                parts = line.strip().split(",")
                if len(parts) < 2:
                    continue
                cluster_id = parts[0]
                for gene_id in parts[1:]:
                    if gene_id:
                        gene_to_cluster[gene_id] = cluster_id

    batch: List[GeneAnnotation] = []
    batch_size = 5000
    count = 0

    with open(annot_file, "r") as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue

            strain_name = parts[1]
            gene_id = parts[5]

            strain = strain_map.get(strain_name)
            cluster_id = gene_to_cluster.get(gene_id)
            cluster = cluster_map.get(cluster_id) if cluster_id else None
            if not strain or not cluster:
                continue

            location = parts[3]
            length = int(parts[4]) if parts[4].isdigit() else None
            gene_name = parts[6] if parts[6] else None
            product = parts[7] if parts[7] else None

            start, end, strand = None, None, None
            loc_match = re.match(r"\[(\d+):(\d+)\]\(([+-])\)", location)
            if loc_match:
                start = int(loc_match.group(1))
                end = int(loc_match.group(2))
                strand = loc_match.group(3)

            batch.append(GeneAnnotation(
                gene_cluster_id=cluster.id,
                strain_id=strain.id,
                gene_id=gene_id,
                product=product,
                gene_name=gene_name,
                contig=parts[2],
                start=start,
                end=end,
                strand=strand,
                length=length,
            ))
            count += 1

            if len(batch) >= batch_size:
                db.add_all(batch)
                db.flush()
                batch = []

    if batch:
        db.add_all(batch)
        db.flush()
    logger.debug("Imported {} gene annotations", count)


# ═══════════════════════════════════════════════════════
# Frequency histogram
# ═══════════════════════════════════════════════════════

def _compute_freq_histogram(
    db: Session, project_id: int,
    cluster_map: Dict[str, GeneCluster],
):
    """Compute cluster-strain frequency histogram from cluster data."""
    freq_counter: Counter = Counter()
    for cluster in cluster_map.values():
        freq_counter[cluster.num_strains] += 1

    batch = [
        ClusterStrainFreq(
            project_id=project_id, num_strains=ns, num_clusters=nc)
        for ns, nc in sorted(freq_counter.items())
    ]
    if batch:
        db.add_all(batch)
        db.flush()
    logger.debug("Computed frequency histogram: {} bins", len(batch))


# ═══════════════════════════════════════════════════════
# Post-processing results
# ═══════════════════════════════════════════════════════

def _import_post_results(db: Session, project_id: int, data_path: Path):
    """Import post stat results (rarefaction, curve_fit, freq, etc.)."""
    post_dir = data_path / "postprocess"
    if not post_dir.is_dir():
        post_dir = data_path

    # Rarefaction
    rarefaction_file = post_dir / "postprocess.rarefaction.tsv"
    if rarefaction_file.exists():
        batch = []
        with open(rarefaction_file, "r") as f:
            f.readline()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                try:
                    batch.append(Rarefaction(
                        project_id=project_id,
                        data_type=parts[0],
                        strain_num=int(parts[1]),
                        value=float(parts[2]),
                    ))
                except (ValueError, IndexError):
                    continue
        if batch:
            db.add_all(batch)
            db.flush()

    # New clusters
    new_clusters_file = post_dir / "postprocess.new_clusters.tsv"
    if new_clusters_file.exists():
        batch = []
        with open(new_clusters_file, "r") as f:
            f.readline()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                try:
                    batch.append(Rarefaction(
                        project_id=project_id,
                        data_type="New clusters",
                        strain_num=int(parts[0]),
                        value=float(parts[1]),
                    ))
                except (ValueError, IndexError):
                    continue
        if batch:
            db.add_all(batch)
            db.flush()

    # Cluster-strain frequency
    freq_file = post_dir / "postprocess.clust_strain_freq.tsv"
    if freq_file.exists():
        freq_counter: Counter = Counter()
        with open(freq_file, "r") as f:
            f.readline()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    try:
                        strain_num = int(parts[1])
                        freq_counter[strain_num] += 1
                    except (ValueError, IndexError):
                        continue
        batch = [
            ClusterStrainFreq(
                project_id=project_id, num_strains=ns, num_clusters=nc)
            for ns, nc in sorted(freq_counter.items())
        ]
        if batch:
            db.add_all(batch)
            db.flush()

    # Curve fit
    curve_file = post_dir / "postprocess.curve_fit.txt"
    if curve_file.exists():
        content = curve_file.read_text().strip()
        if content:
            openness = None
            if "open" in content.lower():
                openness = "open"
            elif "closed" in content.lower():
                openness = "closed"
            db.add(CurveFit(
                project_id=project_id,
                model_name="heaps",
                params_json=content,
                openness=openness,
            ))
            db.flush()

    # Pan group stats
    pan_group_file = post_dir / "postprocess.pan_group_stat.tsv"
    if pan_group_file.exists():
        _import_pan_group_stat(db, project_id, pan_group_file)

    # Paralog stats
    para_stat_file = post_dir / "postprocess.para_stat.tsv"
    if para_stat_file.exists():
        _import_paralog_stat(db, project_id, para_stat_file)

    # Stat attrs
    stat_attrs_file = post_dir / "postprocess.stat_attrs.tsv"
    if stat_attrs_file.exists():
        batch = []
        with open(stat_attrs_file, "r") as f:
            f.readline()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                try:
                    batch.append(StatAttrBin(
                        project_id=project_id,
                        attr=parts[0],
                        group=parts[1],
                        edge=parts[2],
                        count=int(parts[3]),
                    ))
                except (ValueError, IndexError):
                    continue
        if batch:
            db.add_all(batch)
            db.flush()


def _import_pan_group_stat(
    db: Session, project_id: int, pan_group_file: Path,
):
    """Import postprocess.pan_group_stat.tsv → PanGroupStat."""
    group_name_map = {
        "Strict core": "strict_core",
        "Core": "core",
        "Soft core": "soft_core",
        "Shell": "shell",
        "Cloud": "cloud",
    }
    group_stats: Dict[str, Tuple[int, float]] = {}

    with open(pan_group_file, "r") as f:
        f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            category = parts[0].strip()
            try:
                count = int(parts[1].strip())
                pct = float(parts[2].strip())
            except (ValueError, IndexError):
                continue
            group = group_name_map.get(category)
            if group:
                if group in group_stats:
                    old_count, old_pct = group_stats[group]
                    group_stats[group] = (old_count + count, old_pct + pct)
                else:
                    group_stats[group] = (count, pct)

    for group, (count, pct) in group_stats.items():
        db.add(PanGroupStat(
            project_id=project_id,
            group=group,
            num_clusters=count,
            num_genes=count,
            percentage=pct,
        ))
    db.flush()
    logger.debug("Imported pan group stats: {} groups", len(group_stats))


def _import_paralog_stat(
    db: Session, project_id: int, para_stat_file: Path,
):
    """Import postprocess.para_stat.tsv → ParalogStat."""
    batch: List[ParalogStat] = []
    with open(para_stat_file, "r") as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            try:
                batch.append(ParalogStat(
                    project_id=project_id,
                    pan_group=parts[0],
                    para_strain=int(parts[1]),
                    para_gene=int(parts[2]),
                ))
            except (ValueError, IndexError):
                continue
    if batch:
        db.add_all(batch)
        db.flush()
    logger.debug("Imported {} paralog stat rows", len(batch))


# ═══════════════════════════════════════════════════════
# Phylogeny trees
# ═══════════════════════════════════════════════════════

def _import_trees(db: Session, project_id: int, data_path: Path):
    """Import Newick tree files."""
    tree_patterns = [
        ("raw", "*.nwk"),
        ("raw", "*.treefile"),
        ("raw", "*.tree"),
        ("recombination_masked", "*masked*.nwk"),
        ("recombination_masked", "*masked*.treefile"),
    ]
    imported_files: set = set()

    for tree_type, pattern in tree_patterns:
        for tree_file in data_path.rglob(pattern):
            if tree_file in imported_files:
                continue
            if tree_file.is_symlink() and not tree_file.exists():
                continue
            imported_files.add(tree_file)
            try:
                newick = tree_file.read_text().strip()
            except (OSError, IOError):
                continue
            if not newick:
                continue

            actual_type = tree_type
            if ("mask" in tree_file.name.lower()
                    or "gubbins" in tree_file.name.lower()):
                actual_type = "recombination_masked"

            db.add(PhylogenyTree(
                project_id=project_id,
                tree_type=actual_type,
                newick=newick,
                method=_guess_tree_method(tree_file),
            ))
    db.flush()


def _guess_tree_method(tree_file: Path) -> str:
    name_lower = tree_file.name.lower()
    if "fasttree" in name_lower:
        return "fasttree"
    if "iqtree" in name_lower or ".treefile" in name_lower:
        return "iqtree"
    if "raxml" in name_lower:
        return "raxml"
    return "unknown"


# ═══════════════════════════════════════════════════════
# BAPS
# ═══════════════════════════════════════════════════════

def _import_baps(
    db: Session, project_id: int, data_path: Path,
    strain_map: Dict[str, Strain],
):
    """Import fastBAPS clustering results."""
    baps_file = data_path / "fastbaps_clusters.csv"
    if not baps_file.exists():
        candidates = list(data_path.rglob("fastbaps_clusters.csv"))
        if not candidates:
            return
        baps_file = candidates[0]

    with open(baps_file, "r") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if not header:
            return
        for row in reader:
            if len(row) < 2:
                continue
            strain = strain_map.get(row[0])
            if not strain:
                continue
            for level_idx, label in enumerate(row[1:], start=1):
                db.add(BAPSCluster(
                    project_id=project_id,
                    strain_id=strain.id,
                    cluster_label=str(label),
                    level=level_idx,
                ))
    db.flush()


# ═══════════════════════════════════════════════════════
# Preprocess statistics
# ═══════════════════════════════════════════════════════

def _import_prep_stats(
    db: Session, project_id: int, data_path: Path,
    strain_map: Dict[str, Strain],
):
    """Import preprocess.stat.tsv → PrepStat + update Strain fields."""
    stat_file = data_path / "preprocess.stat.tsv"
    if not stat_file.exists():
        return

    batch: List[PrepStat] = []
    with open(stat_file, "r") as f:
        f.readline()  # header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            strain_name = parts[0]
            contig_num = _safe_int(parts[1])
            total_gene_num = _safe_int(parts[2])
            gene_incomplete = _safe_int(parts[3])
            half_core = _safe_int(parts[4])
            single_cloud = _safe_int(parts[5])
            falen = _safe_int(parts[6])
            nucl_comp = parts[7] if parts[7] else None
            ani = _safe_float(parts[8])
            is_darb = _safe_int(parts[9])
            is_outlier_ani = _safe_int(parts[10])
            is_outlier_gene = _safe_int(parts[11])

            ps = PrepStat(
                project_id=project_id,
                strain_name=strain_name,
                contig_num=contig_num,
                total_gene_num=total_gene_num,
                gene_incomplete=gene_incomplete,
                half_core=half_core,
                single_cloud=single_cloud,
                genome_size=falen,
                nucleotide_composition=nucl_comp,
                ani=ani,
                is_darb=is_darb,
                is_outlier_ani=is_outlier_ani,
                is_outlier_gene=is_outlier_gene,
            )
            batch.append(ps)

            # Update strain object with genome info
            strain = strain_map.get(strain_name)
            if strain:
                strain.num_genes = total_gene_num
                strain.genome_size = falen
                strain.num_contigs = contig_num
                if nucl_comp:
                    try:
                        atcg = [int(x) for x in nucl_comp.split("|")]
                        if len(atcg) == 4 and sum(atcg) > 0:
                            strain.gc_content = (atcg[2] + atcg[3]) / sum(atcg)
                    except ValueError:
                        pass

    if batch:
        db.add_all(batch)
        db.flush()
    logger.debug("Imported {} preprocess stat rows", len(batch))


def _import_gene_code(db: Session, project_id: int, data_path: Path):
    """Import preprocess.gene_code.csv → GeneCodeUsage."""
    gc_file = data_path / "preprocess.gene_code.csv"
    if not gc_file.exists():
        return

    batch: List[GeneCodeUsage] = []
    with open(gc_file, "r") as f:
        header = f.readline().strip().split(",")
        strain_names = header[1:]
        for line in f:
            parts = line.strip().split(",")
            if len(parts) < 2:
                continue
            gene_code = parts[0]
            for i, sname in enumerate(strain_names):
                count = _safe_int(parts[i + 1]) if i + 1 < len(parts) else 0
                if count and count > 0:
                    batch.append(GeneCodeUsage(
                        project_id=project_id,
                        gene_code=gene_code,
                        strain_name=sname,
                        count=count,
                    ))
    if batch:
        db.add_all(batch)
        db.flush()
    logger.debug("Imported {} gene code usage entries", len(batch))


# ═══════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════

def _find_file(data_path: Path, *candidates: str) -> Optional[Path]:
    """Return the first existing file from candidates (relative to data_path)."""
    for c in candidates:
        p = data_path / c
        if p.exists():
            return p
    return None


def _safe_int(s: str) -> Optional[int]:
    try:
        return int(s)
    except (ValueError, TypeError):
        return None


def _safe_float(s: str) -> Optional[float]:
    try:
        return float(s)
    except (ValueError, TypeError):
        return None


def parse_readme(readme_path: Path) -> Dict[str, str]:
    """Parse a PGAP2 project README.md for metadata.

    Supports both tab-separated and multi-space-separated key-value pairs.
    """
    KNOWN_KEYS = {
        "project name": "project_name",
        "species": "species",
        "description": "description",
    }
    info: Dict[str, str] = {}
    with open(readme_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("```") or not line:
                continue
            if "\t" in line:
                parts = line.split("\t", 1)
                key = parts[0].strip().lower()
                val = parts[1].strip() if len(parts) > 1 else ""
            else:
                matched = False
                for known_key in sorted(KNOWN_KEYS.keys(), key=len, reverse=True):
                    if line.lower().startswith(known_key):
                        rest = line[len(known_key):]
                        val = rest.strip()
                        key = known_key
                        matched = True
                        break
                if not matched:
                    continue
            field = KNOWN_KEYS.get(key)
            if field and val:
                info[field] = val
    return info


# ═══════════════════════════════════════════════════════
# File embedding — pack disk files into project_files table
# ═══════════════════════════════════════════════════════

# Category → (subdirectory relative to data_path, glob pattern, mime_type)
#
# NOTE: protein (01.gene_prot) and CDS (01.gene_cds) FASTA are intentionally
# NOT embedded.  They are only used by the "Download FASTA" dropdown in
# ClusterMSADialog and would significantly bloat the .db file.  For the
# download-only use case, the Web app falls back to reading them from disk.
_FILE_SPECS: List[Tuple[str, str, str, str]] = [
    # MSA trimmed alignments — displayed in the MSA viewer
    ("msa", "postprocess_phylogeny/04.trim_alignment", "*.fa",
     "text/x-fasta"),
]

# Files searched directly from data_path root
_ROOT_FILE_SPECS: List[Tuple[str, str, str]] = [
    # Annotated layout GML
    ("layout_gml", "pgap2.partition.map_annotated_layout.gml",
     "application/gml+xml"),
    ("layout_gml", "partition/pgap2.partition.map_annotated_layout.gml",
     "application/gml+xml"),
    # Reference-based graph layout JSON (Cytoscape.js)
    ("graph_json", "pgap2.graph_layout.json", "application/json"),
    ("graph_json", "partition/pgap2.graph_layout.json", "application/json"),
    # D3 Canvas precomputed layout JSON
    ("layout_json", "pgap2.partition.map_annotated_layout.json",
     "application/json"),
    ("layout_json", "output/pgap2.partition.map_annotated_layout.json",
     "application/json"),
]


def _embed_files(db: Session, project_id: int, data_path: Path):
    """Walk PGAP2 output dirs and store every Web-needed file as a
    zlib-compressed BLOB in the ``project_files`` table.
    """
    count = 0
    total_raw = 0
    total_compressed = 0
    batch: List[ProjectFile] = []
    batch_size = 50  # flush every N files to limit memory
    category_counts: Dict[str, int] = {}

    # ── Glob-based directories (MSA) ──
    for category, subdir, pattern, mime in _FILE_SPECS:
        full_dir = data_path / subdir
        if not full_dir.is_dir():
            continue
        cat_count = 0
        for fpath in sorted(full_dir.glob(pattern)):
            if not fpath.is_file():
                continue
            rel = fpath.relative_to(data_path).as_posix()
            raw = fpath.read_bytes()
            compressed = zlib.compress(raw, level=6)
            total_raw += len(raw)
            total_compressed += len(compressed)
            batch.append(ProjectFile(
                project_id=project_id,
                path=rel,
                category=category,
                mime_type=mime,
                size=len(raw),
                compressed=True,
                content=compressed,
            ))
            count += 1
            cat_count += 1
            if len(batch) >= batch_size:
                db.add_all(batch)
                db.flush()
                batch = []
        if cat_count:
            category_counts[category] = cat_count

    # ── Root-level files (GML, JSON, TSV) ──
    seen_categories: Dict[str, bool] = {}
    for category, rel_path, mime in _ROOT_FILE_SPECS:
        # For categories with multiple candidates, take the first found
        if category in seen_categories:
            continue
        fpath = data_path / rel_path
        if not fpath.is_file():
            continue
        seen_categories[category] = True
        raw = fpath.read_bytes()
        compressed = zlib.compress(raw, level=6)
        total_raw += len(raw)
        total_compressed += len(compressed)
        batch.append(ProjectFile(
            project_id=project_id,
            path=rel_path,
            category=category,
            mime_type=mime,
            size=len(raw),
            compressed=True,
            content=compressed,
        ))
        count += 1
        category_counts[category] = category_counts.get(category, 0) + 1

    if batch:
        db.add_all(batch)
        db.flush()

    if count:
        ratio = total_raw / total_compressed if total_compressed else 0
        parts = [f"{cat}: {n}" for cat, n in sorted(category_counts.items())]
        logger.info(
            "       {} files embedded ({:.1f} MB -> {:.1f} MB, {:.1f}x compression)",
            count, total_raw / 1024 / 1024,
            total_compressed / 1024 / 1024, ratio,
        )
        logger.info("       [{}]", ", ".join(parts))
    else:
        logger.info("       No files found to embed")
