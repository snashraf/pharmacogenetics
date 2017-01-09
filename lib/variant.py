#!/usr/bin/python
# -*- coding: utf-8 -*-

import json
import urllib2
from modules.pgkb_functions import getRef

# -------------------------------------------------------------------------

class Variant(object):

	"""
	This class stores information on variants. It can take info from pharmkgb,
	 entrez or snpedia (for additional information))
	"""

	def __init__(self, rs):

		# initialize

		self.rs = rs  # rs number of variant

		self.Load()  # load data

		self.GetLocation()

		self.GetAlias()  # get HGVS alias'

		self.SetDefaults()


	def Load(self):

		uri = \
			'https://api.pharmgkb.org/v1/data/variant/?symbol=%s&view=max' \
			% self.rs

		# get data and read in this json file

		data = urllib2.urlopen(uri)

		self.json = json.load(data)[0]

		self.id = self.json['id']

		self.muttype = self.json['type']


	def GetLocation(self):

		if 'GRCh37' not in self.json['location']['name']:

			return

		else:

			self.chr = self.json['location']['name'].split(']'
					)[1].split(':')[0].strip('chr')

			self.begin = self.json['location']['begin']

			self.end = self.json['location']['end']

			self.ref = self.json['location']['reference']

			self.alt = ','.join(self.json['location']['variants'])


	def GetAlias(self):

		# get HGVS alias' for variant, depending on mode again (use later)

		try:

			self.names = self.json['altNames']['synonym']

		except:

			self.names = []

			for doc in self.json['alternateLocations']:

				if 'RefSeq DNA' in doc['sequence']['resource']:

					xref = doc['sequence']['xrefId']

					pos = doc['begin']

					ref = doc['reference']

					alt = doc['variants'][0]

					name = xref + ':g.' + str(pos) + ref + '>' + alt

					self.names.append(name)


	def SetDefaults(self):

		self.nend = self.end
		
		self.nbegin = self.begin
		
		self.nref = self.ref
		
		self.nalt = self.alt

		if self.muttype == "in-del":

			self.LeftShift()


	def LeftShift(self):

		alts = []
		
		# left shift position by 1

		self.nbegin = self.begin - 1

		self.nend = self.end

		# get reference nucleotide at that position
				 
		prevbase = getRef(self.chr, self.nbegin, self.nbegin)

		# set defaults

		if self.ref == "-":

		# scenario 1: insertion (REF - ALT A)

			self.nref = prevbase

			for alt in self.alt.split(","):
				
				alt = prevbase + alt

				alts.append(alt)				 

			self.nalt = ", ".join(alts)

		# scenario 2: deletion ( REF A, ALT -, A)

		elif "-" in self.alt:

			self.nref = prevbase + self.ref

			for alt in self.alt.split(","):

				if alt == "-":

					alt = prevbase

				if alt == self.ref:

					alt = self.nref

				else:

					alts.append(alt)
			
		elif "(" in self.ref:

			# TA repeats

			# manual for now

			self.nref = prevbase + "TA"

			# subtract ref TAs 

			alts.append(prevbase)

			alts.append(prevbase + "TATA")

			alts.append(prevbase + "TATATA")

		self.nalt = ",".join(alts)