"""Preprocess statistics models."""

from typing import Optional

from sqlalchemy import String, Float, Integer, ForeignKey
from sqlalchemy.orm import Mapped, mapped_column

from pgap2.models.base import Base


class PrepStat(Base):
    __tablename__ = "prep_stats"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("projects.id"), nullable=False)
    strain_name: Mapped[str] = mapped_column(String(255), nullable=False)
    contig_num: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    total_gene_num: Mapped[Optional[int]] = mapped_column(
        Integer, nullable=True)
    gene_incomplete: Mapped[Optional[int]] = mapped_column(
        Integer, nullable=True)
    half_core: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    single_cloud: Mapped[Optional[int]] = mapped_column(
        Integer, nullable=True)
    avg_gene_length: Mapped[Optional[float]] = mapped_column(
        Float, nullable=True)
    genome_size: Mapped[Optional[int]] = mapped_column(
        Integer, nullable=True)
    nucleotide_composition: Mapped[Optional[str]] = mapped_column(
        String(64), nullable=True,
        comment="A|T|C|G counts")
    gc_content: Mapped[Optional[float]] = mapped_column(
        Float, nullable=True)
    coding_density: Mapped[Optional[float]] = mapped_column(
        Float, nullable=True)
    ani: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    is_darb: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    is_outlier_ani: Mapped[Optional[int]] = mapped_column(
        Integer, nullable=True)
    is_outlier_gene: Mapped[Optional[int]] = mapped_column(
        Integer, nullable=True)


class GeneCodeUsage(Base):
    __tablename__ = "gene_code_usage"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("projects.id"), nullable=False)
    gene_code: Mapped[str] = mapped_column(String(16), nullable=False)
    strain_name: Mapped[str] = mapped_column(String(255), nullable=False)
    count: Mapped[int] = mapped_column(Integer, default=0)
