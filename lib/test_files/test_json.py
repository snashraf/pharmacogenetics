#!/usr/bin/python
# -*- coding: utf-8 -*-

import urllib2
import json
from tqdm import tqdm
from jinja2 import Template

# --- getting all pairs ---

uri = 'https://api.pharmgkb.org/v1/report/selectPairs'

	# get data and read in this json file

with open("templates/pair.jj", "r") as f:
	template = Template(f.readlines())

data = urllib2.urlopen(uri)

guid = 'nan'

options = 'nan'

results = json.load(data)

for doc in tqdm(results):

	print type(doc)

	sql = template.render(json = doc)

	print sql

# -------------------------