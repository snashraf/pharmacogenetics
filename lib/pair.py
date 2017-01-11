#!/usr/bin/python
# -*- coding: utf-8 -*-

from modules.pgkb_functions import PGKB_connect

# ---------------------------------------------------------------------

class Pair(object):
    '''
Fill in here...
    '''

    def __init__(self, gid, symbol, did, authobj):
    
        self.did = did

        self.gid = gid

        self.symbol = symbol

        self.authobj = authobj

        self.Link()

        self.FindOptions()


    def Link(self):

        results = PGKB_connect(self.authobj, 'clinAnno', self.did, self.gid)

        varids = 'nan'

        if results is None:

            self.varids = "nan"

        elif results is not None:

            self.varids = []

            for doc in results:

                rsid = doc['location']['displayName']

                self.varids.append(rsid)

        if type(varids) == list:

            self.varids = ','.join(list(set(varids)))
        

    def FindOptions(self):

            self.options = "nan"
            
            results = PGKB_connect(self.authobj, 'clinGuide', self.did, self.gid)

            if results is not None:

                self.guid = results['guid']

                optionlist = results['options']

                for gene in optionlist['data']:

                    if self.symbol in gene['symbol']:

                        self.options = ','.join(gene['options'])

            elif results is None:

                self.guid = "nan"

                self.options = "nan"
