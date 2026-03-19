"""Project model — a pan-genome analysis project."""

from typing import Optional

from sqlalchemy import String, Text, Integer
from sqlalchemy.orm import Mapped, mapped_column, relationship

from pgap2.models.base import Base


class Project(Base):
    __tablename__ = "projects"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(255), nullable=False)
    species: Mapped[str] = mapped_column(String(255), nullable=False)
    description: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    num_strains: Mapped[int] = mapped_column(Integer, default=0)
    num_clusters: Mapped[int] = mapped_column(Integer, default=0)
    data_dir: Mapped[str] = mapped_column(String(1024), nullable=False)
    strain_order: Mapped[Optional[str]] = mapped_column(
        Text, nullable=True,
        comment="Comma-separated strain names for PAV column order")

    # Relationships
    strains = relationship(
        "Strain", back_populates="project", cascade="all, delete-orphan")
    gene_clusters = relationship(
        "GeneCluster", back_populates="project", cascade="all, delete-orphan")

    def __repr__(self) -> str:
        return f"<Project(id={self.id}, name={self.name!r})>"
