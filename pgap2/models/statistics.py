"""Statistics-related models (rarefaction, curve fits, pan-group, etc.)."""

from typing import Optional

from sqlalchemy import String, Text, Float, Integer, ForeignKey
from sqlalchemy.orm import Mapped, mapped_column

from pgap2.models.base import Base


class Rarefaction(Base):
    __tablename__ = "rarefaction"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("projects.id"), nullable=False)
    data_type: Mapped[str] = mapped_column(
        String(32), nullable=False,
        comment="e.g. pan, core, new, unique")
    strain_num: Mapped[int] = mapped_column(Integer, nullable=False)
    value: Mapped[float] = mapped_column(Float, nullable=False)


class PanGroupStat(Base):
    __tablename__ = "pan_group_stats"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("projects.id"), nullable=False)
    group: Mapped[str] = mapped_column(String(32), nullable=False)
    num_clusters: Mapped[int] = mapped_column(Integer, default=0)
    num_genes: Mapped[int] = mapped_column(Integer, default=0)
    percentage: Mapped[Optional[float]] = mapped_column(
        Float, nullable=True)


class CurveFit(Base):
    __tablename__ = "curve_fits"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("projects.id"), nullable=False)
    model_name: Mapped[str] = mapped_column(String(64), nullable=False)
    params_json: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    r_squared: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    openness: Mapped[Optional[str]] = mapped_column(
        String(16), nullable=True,
        comment="open or closed")


class ClusterStrainFreq(Base):
    __tablename__ = "cluster_strain_freq"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("projects.id"), nullable=False)
    num_strains: Mapped[int] = mapped_column(Integer, nullable=False)
    num_clusters: Mapped[int] = mapped_column(Integer, nullable=False)


class ParalogStat(Base):
    __tablename__ = "paralog_stats"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("projects.id"), nullable=False)
    pan_group: Mapped[str] = mapped_column(String(32), nullable=False)
    para_strain: Mapped[int] = mapped_column(Integer, default=0)
    para_gene: Mapped[int] = mapped_column(Integer, default=0)


class StatAttrBin(Base):
    __tablename__ = "stat_attr_bins"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("projects.id"), nullable=False)
    attr: Mapped[str] = mapped_column(
        String(64), nullable=False,
        comment="e.g. gene_length, gc_content")
    group: Mapped[str] = mapped_column(String(32), nullable=False)
    edge: Mapped[str] = mapped_column(String(64), nullable=False)
    count: Mapped[int] = mapped_column(Integer, default=0)
