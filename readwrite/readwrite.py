"""
Functions for reading files into GMT and output a GMT object into files
"""
import pygmt as gmt
import operator

def read_gmt(fn, fuzzy=False):
	"""read a txt file into GMT object"""
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
	g = gmt.GMT(d)
	return g

def write_gmt(g, outfn, fuzzy=False):
	d = g.terms
	obj_fuzzy = g.fuzzy
	if obj_fuzzy != fuzzy:
		print 'the output format is in different type of GMT!'
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
	return			

