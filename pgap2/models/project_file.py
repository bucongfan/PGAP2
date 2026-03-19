"""ProjectFile model — embedded binary files stored as compressed BLOBs.

Every file that PGAP2-Web needs at runtime (MSA, FASTA, GML, layout JSON …)
is stored as a zlib-compressed BLOB inside the SQLite database, so the
single ``.db`` file is fully self-contained.
"""

from typing import Optional

from sqlalchemy import String, Text, Integer, LargeBinary, Boolean, ForeignKey
from sqlalchemy.orm import Mapped, mapped_column, relationship

from pgap2.models.base import Base


class ProjectFile(Base):
    __tablename__ = "project_files"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("projects.id"), nullable=False, index=True)

    # Relative path inside the PGAP2 output directory, using '/' separator.
    # Examples:
    #   "postprocess_phylogeny/04.trim_alignment/clust_0001.aln.codon.trimmed.fa"
    #   "pgap2.partition.map.gml"
    #   "pgap2.graph_layout.json"
    path: Mapped[str] = mapped_column(
        String(1024), nullable=False, index=True)

    # Logical category for efficient querying:
    #   "msa", "protein", "cds", "gml", "layout_json", "graph_json",
    #   "annot_tsv", "tree"
    category: Mapped[str] = mapped_column(
        String(64), nullable=False, index=True)

    mime_type: Mapped[Optional[str]] = mapped_column(
        String(128), nullable=True, default="application/octet-stream")

    # Original (uncompressed) size in bytes
    size: Mapped[int] = mapped_column(Integer, nullable=False, default=0)

    # Whether `content` is zlib-compressed
    compressed: Mapped[bool] = mapped_column(
        Boolean, nullable=False, default=True)

    # The file content (zlib-compressed by default)
    content: Mapped[bytes] = mapped_column(LargeBinary, nullable=False)

    project = relationship("Project")

    def __repr__(self) -> str:
        return (f"<ProjectFile(path={self.path!r}, "
                f"category={self.category!r}, size={self.size})>")
