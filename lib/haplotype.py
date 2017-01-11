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

    def __init__(self, gid):

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

        return response

        # TODO catch error

        for doc in response:

            # get various attributes from resulting json file

            starname = doc['name']

            try:

                hgvs = doc['hgvs']

            except KeyError:

                hgvs = None

            # WHY NOT JUST JSON??
            # CONSIDER GENERAL JSON INTErFACE [Get Haplotype etc]
            # CONSIDER JINJA2 TEMPLAETING FOR SQL INTERFACE
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

                    # TODO STOP DUPLICATING DATA

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
