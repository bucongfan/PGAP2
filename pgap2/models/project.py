"""Project model — a pan-genome analysis project."""

from datetime import datetime
from typing import Optional

from sqlalchemy import Boolean, DateTime, String, Text, Integer
from sqlalchemy.orm import Mapped, mapped_column, relationship

from pgap2.models.base import Base


class Project(Base):
    __tablename__ = "projects"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(255), nullable=False, index=True)
    species: Mapped[str] = mapped_column(
        String(255), nullable=False, index=True)
    description: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    species_taxid: Mapped[Optional[int]] = mapped_column(
        Integer, nullable=True)
    num_strains: Mapped[int] = mapped_column(Integer, default=0)
    num_clusters: Mapped[int] = mapped_column(Integer, default=0)
    data_dir: Mapped[str] = mapped_column(String(1024), nullable=False)
    sqlite_path: Mapped[Optional[str]] = mapped_column(
        String(1024), nullable=True)
    strain_order: Mapped[Optional[str]] = mapped_column(
        Text, nullable=True,
        comment="Comma-separated strain names for PAV column order")
    source: Mapped[Optional[str]] = mapped_column(
        String(32), nullable=True, default="database", index=True)
    is_public: Mapped[Optional[bool]] = mapped_column(
        Boolean, nullable=True, default=True)
    owner_id: Mapped[Optional[int]] = mapped_column(
        Integer, nullable=True, index=True)
    job_id: Mapped[Optional[int]] = mapped_column(
        Integer, nullable=True, index=True)
    expires_at: Mapped[Optional[datetime]] = mapped_column(
        DateTime, nullable=True, index=True)
    created_at: Mapped[Optional[datetime]] = mapped_column(
        DateTime, nullable=True, default=datetime.utcnow)
    updated_at: Mapped[Optional[datetime]] = mapped_column(
        DateTime, nullable=True, default=datetime.utcnow,
        onupdate=datetime.utcnow)

    # Relationships
    strains = relationship(
        "Strain", back_populates="project", cascade="all, delete-orphan")
    gene_clusters = relationship(
        "GeneCluster", back_populates="project", cascade="all, delete-orphan")

    def __repr__(self) -> str:
        return f"<Project(id={self.id}, name={self.name!r})>"
