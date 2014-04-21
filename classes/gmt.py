'''
the GMT class
Author: Zichen Wang
3/31/2014

'''
from pygmt import algorithms
from pygmt.readwrite import *

__author__ = "Zichen Wang (wangzc921@gmail.com) "

class GMT(object):

	def __init__(self, data=None, fuzzy=False):
		self.meta = {}
		self.terms = {}
		self.genes = {} # a container for genes with occurrence as value
		self.fuzzy = fuzzy
		if data is not None:
			self.terms, self.fuzzy = to_gmt(data)
			# self.genes = algorithms.count_gene_occ(self.terms)

	@property
	def name(self):
		return self.meta.get('name','')

	@name.setter
	def name(self, s):
		self.meta['name'] = s

	def __str__(self):
		return self.name

	def __iter__(self):
		"""Iterate over the terms. Use the expression 'for t in Gmt:'. 
		Returns an iterator over all terms in the Gmt."""
		return iter(self.terms)

	def __contains__(self, t):
		try:
			return t in self.terms
		except TypeError:
			return False

	def __len__(self):
		return len(self.terms)

	def __getitem__(self, t):
		return self.terms[t]


def to_gmt(data):
	if type(data) == dict: ## data type must be a dictionary
		_key = data.keys()[0]
		if type(data[_key]) != dict:
			fuzzy = False
		else:
			fuzzy = True
		return data, fuzzy
	else:
		raise TypeError("input can't be converted to a GMT")

