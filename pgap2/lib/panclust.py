"""

PanClust is a legacy class that may be deprecated in future versions. 

The primary purpose of this class is to store basic metadata for each gene cluster, including values such as 
minimum, uniqueness, variance, mean, and the representative node. The constructor accepts a set of genes 
and their associated statistics, storing them as attributes of the class.

Additionally, the class maintains a list of per-strain gene clustering records, enabling downstream analyses 
that rely on strain-specific clustering information. PanClust itself does not implement any analysis logic; 
it functions purely as a structured data container.

"""


class Panclust():
    def __init__(self, one_pan: set, min: float, uni: float, var: float, mean: float, repre_node: str) -> None:
        '''
        Record each clutered pan clust

        parameters:
        one_pan: .update({gene: labels[i]})
        uni: uniqueness

        returns:

        '''
        self.min = min
        self.var = var
        self.uni = uni
        self.mean = mean
        self.repre_node = repre_node

        # gene_clust_list = [[] for _ in range(1)]

        # for gene in one_pan:
        #     strain_index = int(one_pan[gene])
        #     gene_clust_list[strain_index].append(gene)
        self.gene_clust_list = one_pan
