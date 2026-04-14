"""ORM models for conserved gene neighborhood analysis."""

from typing import Optional

from sqlalchemy import String, Text, Float, Integer, ForeignKey
from sqlalchemy.orm import Mapped, mapped_column, relationship

from pgap2.models.base import Base


class GeneNeighborhood(Base):
    """A conserved gene neighborhood (syntenic block) identified across genomes."""
    __tablename__ = "gene_neighborhoods"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("projects.id"), nullable=False)
    neighborhood_id: Mapped[str] = mapped_column(
        String(64), nullable=False,
        comment="Display ID, e.g. neigh_0")
    num_clusters: Mapped[int] = mapped_column(
        Integer, default=0,
        comment="Number of distinct gene clusters in this neighborhood")
    num_genomes: Mapped[int] = mapped_column(
        Integer, default=0,
        comment="Number of genomes where this neighborhood is conserved")
    pvalue_cluster: Mapped[Optional[float]] = mapped_column(
        Float, nullable=True,
        comment="Spatial clustering P-value of the representative block")
    pvalue_order: Mapped[Optional[float]] = mapped_column(
        Float, nullable=True,
        comment="Gene-order conservation P-value of the representative block")
    score: Mapped[Optional[float]] = mapped_column(
        Float, nullable=True,
        comment="Combined cluster-match score")
    cluster_ids: Mapped[Optional[str]] = mapped_column(
        Text, nullable=True,
        comment="Comma-separated cluster IDs composing this neighborhood")

    members = relationship(
        "NeighborhoodMember", back_populates="neighborhood",
        cascade="all, delete-orphan")

    def __repr__(self) -> str:
        return (f"<GeneNeighborhood(id={self.id}, "
                f"neighborhood_id={self.neighborhood_id!r}, "
                f"num_clusters={self.num_clusters}, "
                f"num_genomes={self.num_genomes})>")


class NeighborhoodMember(Base):
    """A gene-level member within a conserved gene neighborhood."""
    __tablename__ = "neighborhood_members"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    neighborhood_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("gene_neighborhoods.id"), nullable=False)
    gene_cluster_id: Mapped[Optional[int]] = mapped_column(
        Integer, ForeignKey("gene_clusters.id"), nullable=True)
    strain_id: Mapped[Optional[int]] = mapped_column(
        Integer, ForeignKey("strains.id"), nullable=True)
    contig: Mapped[Optional[str]] = mapped_column(
        String(255), nullable=True)
    position_index: Mapped[Optional[int]] = mapped_column(
        Integer, nullable=True,
        comment="Ordinal gene position on the contig")
    start: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    end: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    strand: Mapped[Optional[str]] = mapped_column(String(1), nullable=True)
    gene_id: Mapped[Optional[str]] = mapped_column(
        String(255), nullable=True,
        comment="Internal gene identifier, e.g. 0:1:23")

    neighborhood = relationship(
        "GeneNeighborhood", back_populates="members")

    def __repr__(self) -> str:
        return (f"<NeighborhoodMember(id={self.id}, "
                f"gene_id={self.gene_id!r})>")
