"""
This module defines the Strain class, which is used to represent a strain in a pangenome analysis.
The Strain class includes attributes such as strain name, strain index, number of genes, and
the number of bad genes. It is primarily used to store metadata about each strain in the context
of pangenome analysis.
It is a simple data structure without any methods for analysis or processing.
"""


class Strain():
    def __init__(self, strain_name, strain_index, bed_gene_num, gene_num) -> None:
        self.strain_name = strain_name
        self.strain_index = strain_index
        self.gene_num = gene_num
        self.bad_gene_num = bed_gene_num
