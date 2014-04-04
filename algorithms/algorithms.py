"""Common algorithms for GMT objects.

Author: Zichen Wang
3/31/2014

"""


def count_gene_occ(d):
	d_gene_occ = {}
	genes = []
	for k in d:
		if type(d[k]) != dict:
			genes += d[k]
		else:
			genes += d[k].keys()
	gene_set = set(genes)
	for gene in gene_set:
		d_gene_occ[gene] = genes.count(gene)
	return d_gene_occ


