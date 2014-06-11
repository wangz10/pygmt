"""
Functions for reading files into GMT and output a GMT object into files
"""
import pygmt as pg
import operator
import warnings

def read_gmt(fn, fuzzy=False, count=False):
	"""read a txt file into GMT object, 
	count specifies whether count the gene frequency"""
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
	g = pg.GMT(d, count=count)
	return g

def write_gmt(g, outfn, fuzzy=False, reverse=True):
	d = g.terms
	obj_fuzzy = g.fuzzy
	if obj_fuzzy != fuzzy:
		warnings.warn(
			'The output GMT format (fuzzy %s) is different from the original type of GMT (fuzzy %s)'%(fuzzy, g.fuzzy),
			RuntimeWarning)
	with open (outfn, 'w') as out:
		for k in d:
			out.write(k + '\tna\t')
			if not fuzzy:
				for gene in d[k]:
					out.write(gene+'\t')
				out.write('\n')
			else:
				assert type(d[k]) == dict
				sorted_d = sorted(d[k].iteritems(), key=operator.itemgetter(1), reverse=reverse)
				for gene, val in sorted_d:
					out.write(gene + ',' + str(val) + '\t')
				out.write('\n')
	return			

def write_df(g, outfn, binarized=True, rorder=None, corder=None, placeholder='0'):
	"""function that write a gmt object into a dataframe format,
	with genes being col names and terms being row names.
	vals are 1/0 if binarized == True; vals are fuzzy values otherwise;
	rorder: (list of terms)
		specifies the subset and the order of rows(terms) in the dataframe
	corder: (list of genes)
		specifies the subset and the order of rows(genes) in the dataframe
		"""
	d = g.terms
	if rorder:
		if len(set(rorder) & set(d)) != len(rorder):
			warnings.warn(
				'Additional terms in rorder do not exist in GMT object, they are filled with 0',
				RuntimeWarning)
	else:
		rorder = d.keys()
	genes = g.genes.keys()
	if corder:
		assert len(set(corder) & set(genes)) == len(corder)
		genes = corder
	with open (outfn, 'w') as out:
		out.write('#\t')
		for gene in genes: ## write col names
			out.write(gene + '\t')
		out.write('\n')
		for term in rorder: ## write rows
			out.write(term + '\t')
			if term in d:
				for gene in genes:
					if binarized:
						if gene in d[term]:
							out.write('1\t')
						else:
							out.write('0\t')
					else:
						if gene in d[term]:
							out.write(str(d[term][gene]) + '\t')
						else:
							out.write(placeholder + '\t')
				out.write('\n')
			else: ## additional terms not in d
				for gene in genes:
					out.write('0\t')
				out.write('\n')
	return



