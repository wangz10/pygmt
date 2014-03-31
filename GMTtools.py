## GMTtools
## created on 3/5/2014 by Zichen Wang

from fileIO import read_gmt, sortD
import operator, random
from scipy.stats import fisher_exact, pearsonr
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.misc import comb

def read_gmt(fn, fuzzy=False):
	d = {}
	with open (fn) as f:
		if not fuzzy:
			for line in f:
				sl = line.strip().split('\t')
				term = sl[0]
				genes = sl[2:]
				d[term] = genes
		else:
			for line in f:
				sl = line.strip().split('\t')
				term = sl[0]
				geneVals = sl[2:]
				d[term] = {}
				for geneVal in geneVals:
					gene, val = geneVal.split(',')
					val = float(val)
					d[term][gene] = val
	return d

def write_gmt(d, outfn, fuzzy=False):
	with open (outfn, 'w') as out:
		for k in d:
			out.write(k + '\tna\t')
			if not fuzzy:
				for gene in d[k]:
					out.write(gene+'\t')
				out.write('\n')
			else:
				assert type(d[k]) == dict
				sorted_d = sorted(d[k].iteritems(), key=operator.itemgetter(1), reverse=True)
				for gene, val in sorted_d:
					out.write(gene + ',' + str(val) + '\t')
				out.write('\n')

def bias_adj(gmt, fuzzy=False):
	# convert a crisp gmt to fuzzy with the 
	# 1/occ as the fuzzy value
	# or original fuzzy val * 1/occ as fuzzy value
	all_genes = []
	for term in gmt:
		genes = list(set(gmt[term]))
		all_genes += genes
	gene_set = set(all_genes)
	d_gene_occ = {}
	for gene in gene_set:
		d_gene_occ[gene] = all_genes.count(gene)
	gmt_fuzzy = {}
	for term in gmt:
		gmt_fuzzy[term] = {}
		for gene in gmt[term]:
			if fuzzy:
				gmt_fuzzy[term][gene] = gmt[term][gene]/float(d_gene_occ[gene])
			else:
				gmt_fuzzy[term][gene] = 1/float(d_gene_occ[gene])
	return gmt_fuzzy

def write_enrichment_matrix(result_dict, dRef, outfn):
	with open (outfn, 'w') as out:
		out.write('\t')
		for term in result_dict: # write header, which are terms from dIn
			out.write(term + '\t')
		out.write('\n')
		for termR in dRef:
			out.write(termR + '\t')
			for term in result_dict:
				out.write(str(result_dict[term][termR]) + '\t')
			out.write('\n')

 
def enrichment_fisher(dRef, dIn, outfn=None, universe=20000):
	result_dict = {}
	for term in dIn:
		result_dict[term] = {}
		genesIn = set(dIn[term])
		b = len(genesIn)
		for termR in dRef:
			genesRef = set(dRef[termR])
			c = len(genesRef)
			a = len(genesRef & genesIn)
			_, p_val = fisher_exact([[a,b],[c,universe]])
			result_dict[term][termR] = p_val
	write_enrichment_matrix(result_dict, dRef, outfn)


def castFuzzy(d_gmt, val=1.0):
	assert type(d_gmt) == dict
	d_gmt_f = {}
	for k in d_gmt:
		if type(d_gmt[k][gene]) == dict:
			return d_gmt
		else:
			d_gmt_f[k] = {}
			for gene in d_gmt[k]:
				d_gmt_f[k][gene] = val
			return d_gmt_f


def gmtT(d_gmt):
	"""transpose gmt file: make genes as terms and 
	terms as genes, only available for fuzzy gmt"""
	d_gmtT = {}
	for term in d_gmt:
		for gene in d_gmt[term]:
			if gene not in d_gmtT:
				d_gmtT[gene] = {}
				d_gmtT[gene][term] = d_gmt[term][gene]
			else:
				d_gmtT[gene][term] = d_gmt[term][gene]
	return d_gmtT


def weight_jaccard(d_gmt, outfn=None):
	"""convert a fuzzy gmt to edge list using fuzzy jaccard distance"""
	d_fs_similarity = {}
	i = 0
	total = comb(len(d_gmt), 2, exact=True)
	for term1, term2 in combinations(d_gmt, 2):
		genes1 = set(d_gmt[term1])
		genes2 = set(d_gmt[term2])
		union = genes1 | genes2
		intersect = genes1 & genes2
		diff1 = genes1 - genes2
		diff2 = genes2 - genes1
		sum_max_weight_intersec = sum([max(d_gmt[term1][gene], d_gmt[term2][gene]) for gene in intersect])
		sum_min_weight_intersec = sum([min(d_gmt[term1][gene], d_gmt[term2][gene]) for gene in intersect])
		sum_weight_diff1 = sum([d_gmt[term1][gene] for gene in diff1])
		sum_weight_diff2 = sum([d_gmt[term2][gene] for gene in diff2])
		jaccard_index = sum_min_weight_intersec/(sum_weight_diff1 + sum_weight_diff2 + sum_max_weight_intersec)
		d_fs_similarity[frozenset([term1,term2])] = jaccard_index
		i += 1
		if i % 10000 ==0:
			print 'processed: ',i,' total: ', total
	print "Finished! Start to sort..."
	sorted_d = sortD(d_fs_similarity)
	print "Finished sorting"
	if outfn:
		with open (outfn, 'w') as out:
			for fs, val in sorted_d:
				for item in fs:
					out.write(item + '\t')
				out.write(str(val) + '\n')
	else:
		return sorted_d

def pearson_corr(d_gmt, outfn=None):
	"""convert a fuzzy gmt to edge list using pearson correlation coefficient"""
	d_fs_similarity = {}
	i = 0
	total = comb(len(d_gmt), 2, exact=True)
	for term1, term2 in combinations(d_gmt, 2):
		genes1 = set(d_gmt[term1])
		genes2 = set(d_gmt[term2])
		assert genes1 == genes2
		vals1 = [d_gmt[term1][gene] for gene in genes1]
		vals2 = [d_gmt[term2][gene] for gene in genes1]
		pcc ,_ = pearsonr(vals1, vals2)
		d_fs_similarity[frozenset([term1,term2])] = pcc
		i += 1
		if i % 10000 ==0:
			print 'processed: ',i,' total: ', total
	print "Finished! Start to sort..."
	sorted_d = sortD(d_fs_similarity)
	print "Finished sorting"
	if outfn:
		with open (outfn, 'w') as out:
			for fs, val in sorted_d:
				for item in fs:
					out.write(item + '\t')
				out.write(str(val) + '\n')
	else:
		return sorted_d





def enrichment_jaccard(dRef, dIn, outfn=None):
	dRef = castFuzzy(dRef)
	dIn = castFuzzy(dIn)
	result_dict = {}
	for term in dIn:
		genesIn = set(dIn[term])
		result_dict[term] = {}
		for termR in dRef:
			genesRef = set(dRef[termR])
			
			intersect = genesIn & genesRef
			union = genesIn | genesRef
			diff1 = genesIn - genesRef
			diff2 = genesRef - genesIn

			sum_min_weight_intersec = 0.0
			sum_max_weight_intersec = 0.0
			sum_weight_diff1 = 0.0
			sum_weight_diff2 = 0.0

			for gene in diff1:
				sum_weight_diff1 += dIn[term][gene]
			for gene in diff2:
				sum_weight_diff2 += dRef[termR][gene]
			for gene in intersec:
				a = dIn[term][gene]
				b = dRef[termR][gene]
				sum_min_weight_intersec += min(a, b)
				sum_max_weight_intersec += max(a, b)
			jaccard = 1 - sum_min_weight_intersec/(sum_weight_diff1 + sum_weight_diff2 + sum_max_weight_intersec)
			result_dict[term][termR] = jaccard
	write_enrichment_matrix(result_dict, dRef, outfn)


def read_enrichment_matrix(filename):
	with open(filename) as f:
		x = []
		for line in f:
			x.append(line)
		termsIn = x[0].strip('\n').split('\t')
		termsIn = termsIn[1:-1]
		d_matrix = {}
		d_termR_val = {}
		termsR = []
		for line in x[1:]:
			y = line.strip().split('\t')
			termR = y[0]
			termsR.append(termR)
			d_termR_val[termR] = y[1:]
		for i in range(len(termsIn)):
			termIn = termsIn[i]
			d_matrix[termIn] = {}
			for termR in d_termR_val:
				d_matrix[termIn][termR] = float(d_termR_val[termR][i])
		return d_matrix


def get_matched_pval_ranks(d_matrix):
	d_termIn_rank = {}
	for termIn in d_matrix:
		d = d_matrix[termIn]
		sorted_d = sorted(d.iteritems(),key=operator.itemgetter(1)) # increasing order
		for rank, (termR, p) in enumerate(sorted_d,start=1):
			if termR == termIn:
				d_termIn_rank[termIn] = rank
	return d_termIn_rank


def plt_DRR(r=None, ax=None, color=None, label=None, ax_right=None, ls=None):
	# r is a list of scaled ranks
	r = np.array(r)
	n = len(r)
	phi = sqrt(1./(4*n))
	sorted_r = np.sort(r)
	
	ecdf = sm.distributions.ECDF(r)
	DR = ecdf(sorted_r)
	DRR = DR - sorted_r
	sDRR = DRR/phi
	print sDRR[-1]

	ax.plot(sorted_r, DRR, color=color, label=label, ls=ls, linewidth=3)
	# ax_right.plot(sorted_r, sDRR, color=color, ls=ls, linewidth=3)
	ax.set_xlabel('Scaled ranks of matched kinases, r', fontsize=24)
	ax.set_xlim([0.0, 1.0])
	ax.set_ylabel('D(r)-r', fontsize=24)	
	ax.set_ylim([-0.02, 0.20])
	ax_right.set_ylabel(r'(D(r)-r)/$\psi$', fontsize=24)
	# to align the ylim of the two axes
	y1, y2 = ax.get_ylim()
	ax_right.set_ylim([y1/phi, y2/phi])
	enlarge_tick_fontsize(ax, 18)
	ax_right.tick_params(axis='both', which='major', labelsize=18)


def get_random_ranks(d_matrix):
	'''permutation of kinase label to get the null distributions'''
	rand_ranks = []
	k1 = d_matrix.keys()[0]
	termsR = d_matrix[k1].keys()
	max_rank = float(len(termsR))
	for termIn in d_matrix:
		d = d_matrix[termIn]
		sorted_d = sorted(d.iteritems(),key=operator.itemgetter(1))
		for i in range(1000):
			termRand = random.choice(termsR)
			for rank, (termR, p) in enumerate(sorted_d,start=1):
				if termR == termRand:
					rand_ranks.append(rank)
	return [x/max_rank for x in rand_ranks]


def plt_random_DRR(r=None, ax=None, color=None, label=None, ax_right=None, ls=None):
	# r is a list of scaled ranks
	r = np.array(r)
	sorted_r = np.sort(r)
	ecdf = sm.distributions.ECDF(r)
	DR = ecdf(sorted_r)
	DRR = DR - sorted_r
	ax.plot(sorted_r, DRR, color=color, label=label, ls=ls, linewidth=3)

def DRR_wrapper():
	return


