#!/usr/bin/python
# -*- coding: utf-8 -*-

import requests
from requests_oauthlib import OAuth2Session
import getpass
import urllib
import urllib2
import ast
import re
import time

# ---------------------------------------------------------------------

def Authenticate():
	"""
	This function creates an authenticating object for use in pharmgkb.
	Necessary for accessing clinical annotations and recommendations!
	:return:
	"""

	print 'Authenticating...'

	val = raw_input("Login as Joanna? y/n: ")

	if "y" in val:

		username = "jwolthuis"

		email = "j.c.wolthuis@students.uu.nl"

	elif "n" in val:

		username = raw_input('PGKB user name: ')

		email = raw_input('PGKB e-mail address: ')

	password = getpass.getpass()

	req = {'username': username, 'email': email, 'password': password}

	# uri for authenticating with this file

	url = 'https://api.pharmgkb.org/v1/auth/oauthSignIn'

	# encode request including user info

	data = urllib.urlencode(req)

	# send request to specified url

	req = urllib2.Request(url, data)

	# read response

	response = urllib2.urlopen(req)

	str_response = response.read()

	# convert token to something that
	# can be used for authentication

	token = ast.literal_eval(str_response)

	# create an authorized session and
	# return this session

	client = OAuth2Session(token=token, auto_refresh_url=url)

	return client


def getJson(uri, client):
	"""
	Uses uri and given authenticated client to get json data.
	:param uri: filled in uri refining query
	:param client: result of Authenticate()
	:return:
	"""

	resp = requests.head(uri)

	status = resp.status_code

	if status == 200:

		r = client.get(uri)

		data = r.json()

		return data

	else:

		return None


def seqMaker(rsidorder, reference, rsids):

	seq = ''

	for rsid in rsidorder:

		try:
			base = rsids[rsid]

		except KeyError:

			base = reference[rsid]

		seq += base

	return seq


def getRef(loc, start, end):

	server = 'http://grch37.rest.ensembl.org'

	ext = '/sequence/region/human/{}:{}..{}?'.format(loc.lstrip('chr'),
			start, end)

	r = requests.get(server + ext, headers={'Content-Type': 'text/plain'
					 })

	if not r.ok:

		r.raise_for_status()

		sys.exit()

	time.sleep(0.3)

	return r.text


def hg19conv(rsid, alt, gid):

	d = {}

	d['id'] = rsid

	d['gid'] = gid

	# find chromosome number

	d['loc'] = rsid.split(':')[0]

	# find location

	d['begin'] = int(rsid[rsid.find(':') + 1:rsid.find('(')])

	d['end'] = d['begin']

	# get reference position

	d['ref'] = getRef(d['loc'], d['begin'], d['begin'])

	d['alt'] = alt

	if 'delGENE' in alt:

		return None

	if d['ref'] == alt:

		return None

	if len(d['ref']) == len(alt):

		d['muttype'] = 'snp'

	if 'ins' in alt:

		d['begin'] += 1

		d['ref'] = "-"

		d['alt'] = alt.replace('ins', "")

		d['muttype'] = 'in-del'

		d['end'] = d['begin'] + len(alt)

	elif 'del' in alt:

		d['alt'] = "-"

		d['muttype'] = 'in-del'

		d['end'] = d['begin'] + len(d['ref'])
	return d

def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

def DeconstructGuideline(markdown):
	markdown = markdown.split("\n")
	tex = ""
	# ------------
	conv = {
"%":"\%",
	}
	# ------------
	nontable = []
	table = []
	for line in markdown:
		if "|" in line:
			table.append(line)
		else:
			nontable.append(line)
	if len(table) == 0:
			tex += unicode(line + "\\\\")
	else:
		for i, line in enumerate(table):
			print line
			spl_line = [item for item in line.split("|") if line != ""]
			andjoin = "&".join(spl_line)
			tex += unicode(andjoin + "\\\\")

	return tex
