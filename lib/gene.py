#!/usr/bin/python
# -*- coding: utf-8 -*-

import urllib2
import json

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

        uri = 'https://api.pharmgkb.org/v1/data/gene/{}?view=max'.format(self.gid)

        try:

            response = urllib2.urlopen(uri)

        except urllib2.HTTPError:
            # TODO Log error message to stderr
            return None

        self.json = json.load(response)