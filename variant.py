#!/usr/bin/python
# -*- coding: utf-8 -*-

import json
import urllib2


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

    def Load(self):

        # check for which mode is being used and adjust query accordingly

        uri = \
            'https://api.pharmgkb.org/v1/data/variant/?symbol=%s&view=max' \
            % self.rs

        # get data and read in this json file

        data = urllib2.urlopen(uri)

        self.json = json.load(data)[0]

        self.id = self.json['id']

        self.type = self.json['type']

        return


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

        if self.mode == 'pharmgkb':

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

class Indel(Variant):
	'''
lalalalalala
	'''

	def __init__(self, rsid):

		pass

