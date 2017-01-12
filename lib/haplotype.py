#!/usr/bin/python
# -*- coding: utf-8 -*-

import urllib2
import json
import modules.pgkb_functions


# -------------------------------------------------------------------------

class Gene:

    """
This class gets information on a given gene based on gene ID. Can use either PharmGKB or Entrez.
    """

    def __init__(self, hapid):

        # initialize with geneID

        self.gid = gid

        self.haplotypes = []

        self.genes = []

        self.Load()  # loads data on gene


    def Load(self):

        # get a list of known haplotype IDs from gene ID

        uri = 'https://api.pharmgkb.org/v1/data/haplotype?gene.accessionId={}&view=max'.format(self.gid)

        data = urllib2.urlopen(uri)

        response = json.load(data)

        # TODO catch error

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

            for var in alleles:

                try:

                    rsid = var['location']['displayName']

                    alt = var['allele']

                    rsids.append((rsid, alt))
