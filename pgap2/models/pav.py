"""Presence-absence variation model."""

from sqlalchemy import Integer, ForeignKey
from sqlalchemy.orm import Mapped, mapped_column, relationship

from pgap2.models.base import Base


class PAV(Base):
    __tablename__ = "pav"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("projects.id"), nullable=False)
    gene_cluster_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("gene_clusters.id"), nullable=False)
    strain_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("strains.id"), nullable=False)
    copy_number: Mapped[int] = mapped_column(Integer, default=1)

    gene_cluster = relationship("GeneCluster", back_populates="pav_entries")
    strain = relationship("Strain")

    def __repr__(self) -> str:
        return (f"<PAV(cluster={self.gene_cluster_id}, "
                f"strain={self.strain_id}, cn={self.copy_number})>")
