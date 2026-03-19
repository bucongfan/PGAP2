"""Gene-cluster model and PanGroup enum."""

import enum
from typing import Optional

from sqlalchemy import String, Text, Float, Integer, Enum, ForeignKey
from sqlalchemy.orm import Mapped, mapped_column, relationship

from pgap2.models.base import Base


class PanGroup(str, enum.Enum):
    """Pan-genome group classification."""
    STRICT_CORE = "strict_core"
    CORE = "core"
    SOFT_CORE = "soft_core"
    SHELL = "shell"
    CLOUD = "cloud"


class GeneCluster(Base):
    __tablename__ = "gene_clusters"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("projects.id"), nullable=False)
    cluster_id: Mapped[str] = mapped_column(String(64), nullable=False)
    annotation: Mapped[Optional[str]] = mapped_column(
        String(1024), nullable=True)
    num_strains: Mapped[int] = mapped_column(Integer, default=0)
    num_genes: Mapped[int] = mapped_column(Integer, default=0)
    group: Mapped[Optional[PanGroup]] = mapped_column(
        Enum(PanGroup), nullable=True)
    frequency: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    avg_length: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    pav_vector: Mapped[Optional[str]] = mapped_column(
        Text, nullable=True,
        comment="Comma-separated 0/1/n copy-number string")

    project = relationship("Project", back_populates="gene_clusters")
    pav_entries = relationship(
        "PAV", back_populates="gene_cluster", cascade="all, delete-orphan")
    gene_annotations = relationship(
        "GeneAnnotation", back_populates="gene_cluster",
        cascade="all, delete-orphan")

    def __repr__(self) -> str:
        return f"<GeneCluster(id={self.id}, cluster={self.cluster_id!r})>"
