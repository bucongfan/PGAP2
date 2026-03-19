"""Phylogeny tree model."""

from typing import Optional

from sqlalchemy import String, Text, Integer, ForeignKey
from sqlalchemy.orm import Mapped, mapped_column

from pgap2.models.base import Base


class PhylogenyTree(Base):
    __tablename__ = "phylogeny_trees"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    project_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("projects.id"), nullable=False)
    tree_type: Mapped[str] = mapped_column(
        String(32), nullable=False, index=True,
        comment="e.g. single_copy, pav")
    newick: Mapped[str] = mapped_column(Text, nullable=False)
    num_leaves: Mapped[int] = mapped_column(Integer, default=0)
    method: Mapped[Optional[str]] = mapped_column(String(64), nullable=True)

    def __repr__(self) -> str:
        return (f"<PhylogenyTree(id={self.id}, "
                f"type={self.tree_type!r})>")
