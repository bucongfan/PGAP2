"""Tajima's D model."""

from typing import Optional

from sqlalchemy import String, Float, Integer, ForeignKey
from sqlalchemy.orm import Mapped, mapped_column

from pgap2.models.base import Base


class TajimasD(Base):
    __tablename__ = "tajimas_d"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("projects.id"), nullable=False)
    cluster_id: Mapped[str] = mapped_column(String(64), nullable=False)
    num_sequences: Mapped[int] = mapped_column(Integer, default=0)
    segregating_sites: Mapped[int] = mapped_column(Integer, default=0)
    tajimas_d: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    pi: Mapped[Optional[float]] = mapped_column(Float, nullable=True)

    def __repr__(self) -> str:
        return (f"<TajimasD(cluster={self.cluster_id!r}, "
                f"d={self.tajimas_d})>")
