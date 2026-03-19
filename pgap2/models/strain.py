"""Strain model — one genome / assembly in a project."""

from typing import Optional

from sqlalchemy import String, Float, Integer, ForeignKey
from sqlalchemy.orm import Mapped, mapped_column, relationship

from pgap2.models.base import Base


class Strain(Base):
    __tablename__ = "strains"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("projects.id"), nullable=False)
    name: Mapped[str] = mapped_column(String(255), nullable=False)
    assembly_accession: Mapped[Optional[str]] = mapped_column(
        String(64), nullable=True)
    num_genes: Mapped[int] = mapped_column(Integer, default=0)
    genome_size: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    gc_content: Mapped[Optional[float]] = mapped_column(
        Float, nullable=True)
    num_contigs: Mapped[Optional[int]] = mapped_column(
        Integer, nullable=True)
    coding_density: Mapped[Optional[float]] = mapped_column(
        Float, nullable=True)

    project = relationship("Project", back_populates="strains")

    def __repr__(self) -> str:
        return f"<Strain(id={self.id}, name={self.name!r})>"
