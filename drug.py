import urllib2
import json

class Drug(object):
	'''
	
	'''
	def __init__(self, did):

	    uri = \
	        'https://api.pharmgkb.org/v1/data/chemical/%s?view=max' \
	        % did

	    self.load()

	def Load(self):

	    # get data and read in this json file

	    data = urllib2.urlopen(uri)

	    json = json.load(data)

	    self.name = json['name']

	    self.terms = json['terms']