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


def cast_fuzzy(gmt, val=1.0):
	""" cast a crisp gmt object to a fuzzy one by
	assigning an arbitury value 1.0 by default.
	"""
	assert type(gmt.fuzzy) == False
	terms_original = gmt.terms
	terms_new = {}
	for term in terms_original:
		terms_new[term] = {}
		for gene in terms_original[k]:
			terms_new[term][gene] = val
	gmt.terms = terms_new
	gmt.fuzzy = True
	return gmt


def gmtT(gmt):
	"""transpose gmt object: make genes as terms and 
	terms as genes, only available for fuzzy gmt"""
	terms_original = gmt.terms
	termsT = {}
	for term in terms_original:
		for gene in terms_original[term]:
			if gene not in termsT:
				termsT[gene] = {}
				termsT[gene][term] = terms_original[term][gene]
			else:
				termsT[gene][term] = terms_original[term][gene]
	gmt.terms = termsT
	return gmt

def union(gmt1, gmt2):
	"""union two gmt objects, may wanna extend to more..."""
	d1 = gmt1.term
	d2 = gmt2.term
	d1.update(d2)
	d_gene_occ = count_gene_occ(d1)

	gmt = GMT(d1)
	gmt.genes = d_gene_occ

	return gmt
	

