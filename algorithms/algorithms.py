"""Common algorithms for GMT objects.

Author: Zichen Wang
3/31/2014

"""
import operator
from collections import Counter
from scipy.stats import fisher_exact
from math import floor
import pygmt as pg

def sortD(d, reverse=True):
	"""sort dict by value in decreasing order (default)"""
	assert type(d) == dict
	sorted_d = sorted(d.iteritems(), key=operator.itemgetter(1), reverse=reverse)
	return sorted_d

def count_gene_occ(d):
	"""return a dictionary of gene and their occurrences"""
	d_gene_occ = {}
	genes = []
	for k in d:
		if type(d[k]) != dict:
			genes += d[k]
		else:
			genes += d[k].keys()
	c = Counter(genes)
	d_gene_occ = dict(c)
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
	d1 = gmt1.terms
	d2 = gmt2.terms
	d1.update(d2)
	d_gene_occ = count_gene_occ(d1)

	gmt = GMT(d1)
	gmt.genes = d_gene_occ

	return gmt
	
def subset(gmt, terms):
	"""subset a gmt objects by a list/set of terms"""
	d = gmt.terms
	terms_sub = dict((term, d[term]) for term in terms)
	gmt.terms = terms_sub
	return gmt

def unify_length(gmt, length=None, reverse=True):
	"""keep top genes for each gmt lines in decreasing order (default)
	of the fuzzy value """
	d = gmt.terms
	terms_new = {}
	for term in d:
		sorted_d = sorted(d[term].iteritems(), key=operator.itemgetter(1), reverse=reverse)
		if length <= len(sorted_d):
			terms_new[term] = dict(sorted_d[0:length])
		else:
			terms_new[term] = dict(sorted_d)
	g = pg.GMT()
	g.terms = terms_new
	g.genes = count_gene_occ(terms_new)
	return g

def remove_by_occ(gmt, percent):
	"""remove genes in the gmt if the occurrences of genes lie in the
	percentile specified:
	e.g.: percent = 10 # remove top 10% """
	d_gene_occ = gmt.genes
	sorted_gene_occ = sortD(d_gene_occ)
	num_rm_genes = floor(len(d_gene_occ) * percent * 0.01)
	genes_to_rm = []
	for i, (gene, occ) in enumerate(sorted_gene_occ, start=1):
		if i != num_rm_genes:
			genes_to_rm.append(gene)
		else:	
			break
	d = gmt.terms
	d_new = {}
	for term in d:
		d_new[term] = [gene for gene in d[term] if gene not in genes_to_rm]
		if len(d_new[term]) == 0:
			del d_new[term]
	g = pg.GMT()
	g.terms = d_new
	return g


def cross_enrichment(gmt1, gmt2, method='fisher',outfn=None, universe=20000):
	"""a function performing enrichment analysis for pairwise terms in two gmts"""
	with open (outfn, 'w') as out:
		out.write('#\t') ## write header
		for term2 in gmt2.terms: 
			out.write(term2+'\t') 
		out.write('\n')
		i = 1
		for term1 in gmt1.terms:
			out.write(term1 + '\t')
			genes1 = set(gmt1.terms[term1])
			print 'processing ', term1, i
			i += 1
			for term2 in gmt2.terms:
				genes2 = set(gmt2.terms[term2])
				a = len(genes1 & genes2)
				b = len(genes1)
				c = len(genes2)
				if method == 'fisher':
					_, p_val = fisher_exact([[a,b],[c,universe]])
				out.write(str(p_val) + '\t')
			out.write('\n')
	return

