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

    def __init__(self, geneid):

        # initialize with geneID

        self.gid = geneid

        self.alleles = []

        self.genes = []

        self.Load()  # loads data on gene

        # self.GetDrugs()  # loads info on associated drugs
        # self.GetDesc()  # gets summary of gene

        try:

            self.GetHaps()

        except urllib2.HTTPError:
            
            pass


    def Load(self):

        uri = 'https://api.pharmgkb.org/v1/data/gene/%s?view=max' \
                % self.gid

        print uri

        try:

            response = urllib2.urlopen(uri)

        except urllib2.HTTPError:

            return None

        self.json = response.load(data)

        self.name = self.json['symbol']

        self.chr = self.json['chr']['name'].lstrip('chr')

        self.start = self.json['chrStart']

        self.stop = self.json['chrStop']


    def GetHaps(self):

        # get a list of known haplotype IDs from gene ID

        uri = \
            'https://api.pharmgkb.org/v1/data/haplotype?gene.accessionId=%s&view=max' \
            % self.gid

        data = urllib2.urlopen(uri)

        response = json.load(data)

        for doc in response:

            # get various attributes from resulting json file

            starname = doc['name']

            try:

                hgvs = doc['hgvs']

            except KeyError:

                hgvs = None

            alleles = doc['alleles']

            hapid = doc['id']

            copynum = doc['copyNumber']

            rsids = []

            dic = {}

            # get all the involved rsids

            for allele in alleles:

                try:

                    rsid = allele['location']['displayName']

                    alt = allele['allele']

                    rsids.append((rsid, alt))
                except KeyError:

                    continue

            guideline = False

            # add to alleles dictionary

            all = {
                'starname': starname,
                'hgvs': hgvs,
                'id': hapid,
                'copynum': copynum,
                'rsids': rsids,
                'guideline': guideline,
                }
            
            self.alleles.append(all)
