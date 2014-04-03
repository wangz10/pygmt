'''
the GMT class
Author: Zichen Wang
3/31/2014

'''

__author__ = "Zichen Wang (wangzc921@gmail.com) "

class GMT(object):

	def __init__(self, data=None, fuzzy=False):
		self.meta = {}
		self.terms = {}
		self.fuzzy = fuzzy
		if data is not None:
			self.terms, self.fuzzy = to_gmt(data)

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
	if type(data) == dict:
		for key in data:
			if type(data[key]) != dict:
				fuzzy = False
			else:
				fuzzy = True
		
		return data, fuzzy
	else:
		raise TypeError("input can't be converted to a GMT")


## testing block below:

d = dict(a=1,b=2)
g = GMT(d)
g2 = GMT(1)
# print g
print g.terms
print g.fuzzy