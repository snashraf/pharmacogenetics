#!/usr/bin/python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------


class Drug(object):
	'''

	'''
	def __init__(self, did):

	    self.uri = \
	        'https://api.pharmgkb.org/v1/data/chemical/%s?view=max' \
	        % did

	    self.Load()


	def Load(self):

	    # get data and read in this json file

	    data = urllib2.urlopen(self.uri)

	    self.json = json.load(data)