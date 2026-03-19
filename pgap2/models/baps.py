"""BAPS clustering model."""

from sqlalchemy import String, Integer, ForeignKey
from sqlalchemy.orm import Mapped, mapped_column

from pgap2.models.base import Base


class BAPSCluster(Base):
    __tablename__ = "baps_clusters"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("projects.id"), nullable=False)
    strain_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("strains.id"), nullable=False)
    cluster_label: Mapped[str] = mapped_column(String(64), nullable=False)
    level: Mapped[int] = mapped_column(Integer, default=1)

    def __repr__(self) -> str:
        return (f"<BAPSCluster(strain={self.strain_id}, "
                f"label={self.cluster_label!r})>")
