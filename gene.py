#!/usr/bin/python
# -*- coding: utf-8 -*-

import urllib2
import json
import pgkb_functions


# -------------------------------------------------------------------------

class Gene:

    """
This class gets information on a given gene based on gene ID. Can use either PharmGKB or Entrez.
    """

    def __init__(self, gid):

        # initialize with geneID

        self.gid = gid

        self.haplotypes = []

        self.genes = []

        self.Load()  # loads data on gene

        # self.GetDrugs()  # loads info on associated drugs
        # self.GetDesc()  # gets summary of gene

        try:

            self.GetHaps()

        except urllib2.HTTPError:
            
            pass


    def Load(self):

        uri = 'https://api.pharmgkb.org/v1/data/gene/{}?view=max'.format(self.gid)

        try:

            response = urllib2.urlopen(uri)

        except urllib2.HTTPError:

            return None

        self.json = json.load(response)

        self.name = self.json['symbol']

        self.chr = self.json['chr']['name'].lstrip('chr')

        self.start = self.json['chrStart']

        self.stop = self.json['chrStop']


    def GetHaps(self):

        # get a list of known haplotype IDs from gene ID

        uri = \
            'https://api.pharmgkb.org/v1/data/haplotype?gene.accessionId={}&view=max' \
            .format(self.gid)

        data = urllib2.urlopen(uri)

        response = json.load(data)

        for doc in response:

            # get various attributes from resulting json file

            starname = doc['name']

            try:

                hgvs = doc['hgvs']

            except KeyError:

                hgvs = None

            haplotypes = doc['alleles']

            hapid = doc['id']

            copynum = doc['copyNumber']

            rsids = []

            dic = {}

            # get all the involved rsids

            for hap in haplotypes:

                try:

                    rsid = hap['location']['displayName']

                    alt = hap['allele']

                    rsids.append((rsid, alt))

                except KeyError:

                    continue

            guideline = False

            # add to alleles dictionary

            d = {
                
                'starname': starname,
                
                'hgvs': hgvs,
                
                'id': hapid,
                
                'copynum': copynum,
                
                'rsids': rsids,
                
                'guideline': guideline,

                }
            
            self.haplotypes.append(d)
