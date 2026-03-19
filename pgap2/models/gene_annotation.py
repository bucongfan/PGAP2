"""Gene annotation model."""

from typing import Optional

from sqlalchemy import String, Integer, ForeignKey
from sqlalchemy.orm import Mapped, mapped_column, relationship

from pgap2.models.base import Base


class GeneAnnotation(Base):
    __tablename__ = "gene_annotations"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    gene_cluster_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("gene_clusters.id"), nullable=False)
    strain_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("strains.id"), nullable=False)
    gene_id: Mapped[str] = mapped_column(String(255), nullable=False)
    product: Mapped[Optional[str]] = mapped_column(String(512), nullable=True)
    gene_name: Mapped[Optional[str]] = mapped_column(
        String(64), nullable=True)
    locus_tag: Mapped[Optional[str]] = mapped_column(
        String(64), nullable=True)
    contig: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    start: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    end: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    strand: Mapped[Optional[str]] = mapped_column(String(1), nullable=True)
    length: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)

    gene_cluster = relationship(
        "GeneCluster", back_populates="gene_annotations")
    strain = relationship("Strain")

    def __repr__(self) -> str:
        return f"<GeneAnnotation(gene_id={self.gene_id!r})>"
