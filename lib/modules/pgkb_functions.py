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

	client = OAuth2Session(token=token)

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


def PGKB_connect(
	authobj,
	mode,
	a,
	b,
	):
	"""
	Find connections between two pharmgkb objects.
	Usually a haplotype or rsid compared with a drug id.
	:param authobj: authentication session resulting from Authenticate()
	:param a: object 1 to compare
	:param did: object 2 to compare, usually a drug id
	:return:
	"""

	if mode == 'clinAnno':

		uri = \
			'https://api.pharmgkb.org/v1/report/pair/{}/{}/clinicalAnnotation?view=base' \
			.format(a, b)

		result = getJson(uri, authobj)

		return result
	
	elif mode == 'clinGuide':

		sources = ['cpic', 'dpwg', 'pro']

		for source in sources:

			# uses did and gid

			uri = \
				'https://api.pharmgkb.org/v1/data/guideline?source={}&relatedChemicals.accessionId={}&relatedGenes.accessionId={}&view=max' \
				.format(source, a, b)

			guidelines = getJson(uri, authobj)

			guid = 'nan'

			options = 'nan'

			if guidelines == None:

				continue

			else:

				for doc in guidelines:

					if doc['objCls'] == 'Guideline':

						guid = doc['id']

						uri = \
							'https://api.pharmgkb.org/v1/report/guideline/{}/options' \
							.format(guid)

						doseOptions = getJson(uri, authobj)

						result = {'guid': guid, 'options': doseOptions}

						return result


def seqMaker(rsidorder, reference, rsids):

	seq = ''

	for rsid in rsidorder:

		try:

			base = rsids[rsid]

		except KeyError:

			base = reference[rsid]

		if 'del' in base or base == '-':

			continue

		if ',' in base:

			bases = base.split(',')

			base = bases[0]

		elif '[' in base or '(' in base:

			# find bracketed objects

			filt = re.split('\[(.*?)\]|\((.*?)\)', base)

			for frag in filt:

				if frag is not None:

					try:

						num = int(frag)

					except:

						motif = frag

			base = motif * num

		seq += base
			
	return seq


def Aligner(seqlist):
	
	pass


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

	# find chromosome number

	d['loc'] = rsid.split(':')[0]

	# find location

	d['begin'] = int(rsid[rsid.find(':') + 1:rsid.find('(')])

	d['end'] = d['begin']

	# get reference position

	d['ref'] = getRef(d['loc'], d['begin'], d['begin'])

	d['alt'] = alt

	if d['ref'] == alt or 'delGENE' in alt:

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

	else:

		return None

	print d

	return d