"""Strain model — one genome / assembly in a project."""

from datetime import datetime
from typing import Optional

from sqlalchemy import DateTime, String, Text, Float, Integer, ForeignKey
from sqlalchemy.orm import Mapped, mapped_column, relationship

from pgap2.models.base import Base


class Strain(Base):
    __tablename__ = "strains"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("projects.id"), nullable=False)
    name: Mapped[str] = mapped_column(String(255), nullable=False, index=True)
    assembly_accession: Mapped[Optional[str]] = mapped_column(
        String(64), nullable=True, index=True)
    num_genes: Mapped[Optional[int]] = mapped_column(
        Integer, nullable=True)
    num_cds: Mapped[Optional[int]] = mapped_column(
        Integer, nullable=True)
    total_genes: Mapped[Optional[int]] = mapped_column(
        Integer, nullable=True)
    genome_size: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    gc_content: Mapped[Optional[float]] = mapped_column(
        Float, nullable=True)
    n50: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    avg_gene_length: Mapped[Optional[float]] = mapped_column(
        Float, nullable=True)
    completeness: Mapped[Optional[float]] = mapped_column(
        Float, nullable=True)
    contamination: Mapped[Optional[float]] = mapped_column(
        Float, nullable=True)
    gene_count: Mapped[Optional[int]] = mapped_column(
        Integer, nullable=True)
    num_contigs: Mapped[Optional[int]] = mapped_column(
        Integer, nullable=True)
    gff_path: Mapped[Optional[str]] = mapped_column(
        String(1024), nullable=True)
    extra_meta: Mapped[Optional[str]] = mapped_column(
        Text, nullable=True, comment="JSON blob for additional metadata")
    coding_density: Mapped[Optional[float]] = mapped_column(
        Float, nullable=True)
    annotation_source: Mapped[Optional[str]] = mapped_column(
        String(1024), nullable=True)
    annotation_json: Mapped[Optional[str]] = mapped_column(
        Text, nullable=True)
    created_at: Mapped[Optional[datetime]] = mapped_column(
        DateTime, nullable=True, default=datetime.utcnow)
    updated_at: Mapped[Optional[datetime]] = mapped_column(
        DateTime, nullable=True, default=datetime.utcnow,
        onupdate=datetime.utcnow)

    project = relationship("Project", back_populates="strains")

    def __repr__(self) -> str:
        return f"<Strain(id={self.id}, name={self.name!r})>"
