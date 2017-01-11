#!/usr/bin/python
# -*- coding: utf-8 -*-

import urllib2
import json
from tqdm import tqdm
from jinja2 import Template
import re

# -------------------------

def createSQL(path, context):
	
	conv = ""
	filt = '\t\n'
	
	with open(path, "r") as f:
		for line in f.readlines():
			conv += re.sub('[\n\t]', '', line)
	
	template = Template(conv)

	result = template.render(json = context)

	return result

# --- getting all pairs ---

uri = 'https://api.pharmgkb.org/v1/report/selectPairs'

	# get data and read in this json file

data = urllib2.urlopen(uri)

results = json.load(data)

for doc in tqdm(results):

	sql = createSQL("templates/insert/pair.jj", doc)

	print sql

# -------------------------