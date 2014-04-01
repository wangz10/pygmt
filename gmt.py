'''
the Gmt class
Author: Zichen Wang
3/31/2014

'''
__author__ = "Zichen Wang (wangzc921@gmail.com) "

class GMT(object):

	def __init__(self, data=None, fuzzy=False):
		self.terms = {}
		if data is not None:
			self.terms = convert.to_gmt(data) ## convert.py to be written for functions making GMT
