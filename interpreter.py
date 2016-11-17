#!/usr/bin/python
# -*- coding: utf-8 -*-
from patient import Patient
from data import DataCollector
import json
import sqlite3
import requests
from requests_oauthlib import OAuth2Session
import urllib
import urllib2
import ast


# -------------------------------------------------------------------------

def Authenticate():
    """
    This function creates an authenticating object for use in pharmgkb.
    Necessary for accessing clinical annotations and recommendations!
    :return:
    """
    print 'Authenticating...'
    req = ''

    # load file with authenticating info, use to request a key

    with open('config/auth.txt', 'r') as f:
        req = json.load(f)

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
    r = client.get(uri)
    data = r.json()

    # if there is just a dict and not a list it is an error message

    if type(data) == dict:
        return None

    # otherwise return valid results

    elif type(data) == list:
        return data


def PGKB_connect(authobj, a, did):
    """
    Find connections between two pharmgkb objects.
    Usually a haplotype or rsid compared with a drug id.
    :param authobj: authentication session resulting from Authenticate()
    :param a: object 1 to compare
    :param did: object 2 to compare, usually a drug id
    :return:
    """
    uri = \
        'https://api.pharmgkb.org/v1/report/pair/%s/%s/clinicalAnnotation' \
        % (a, did)

    # get json file from filled in uri, use authentication client

    result = getJson(uri, authobj)
    results = []

    # collect matching results

    if result is not None:
        for doc in result:
            for phen in doc['allelePhenotypes']:
                results.append(phen)
    else:
        return
    return results

# -------------------------------------------------------------------------

class Interpreter:


    def __init__(self, patientObj):
        self.p = patientObj
        self.advice = {}
        self.conn = sqlite3.connect('pharmacogenetics.db')
        self.sql = self.conn.cursor()

    def DrugAdvice(self):

        authobj = Authenticate()

        # from haplotypes = [PA..., PA...]
        # from rsids = [rs..., rs..., rs...]

        self.sql.execute('''SELECT d.rsid, p.chr1, p.chr2
                            FROM patientvars p
                            JOIN design d ON d.pos = p.pos
                            JOIN variants v ON d.rsid = v.rsid''')

        var_alleles = {rsid:( chr1 + chr2 )
        for( rsid, chr1, chr2 )
        in self.sql.fetchall() }

        input = self.p.haplotypes + var_alleles.keys()

        for var in input:

            # check for haplotype or rs nomenclature

            if "PA" in var:
                self.sql.execute("SELECT DISTINCT gid FROM alleles WHERE hapid = ?", (var,))
                var_alleles[var] = self.p.starhaps[var]

            elif "rs" in var:
                self.sql.execute("SELECT DISTINCT gid FROM variants WHERE rsid = ?", (var,))

            gids = [tup[0] for tup in self.sql.fetchall()]

            # find gene ids associated with this variant or haplotype

            for gid in gids:

                # for each gene, find associated drugs

                self.sql.execute("SELECT DISTINCT did FROM drugpairs WHERE gid = ?", (gid,))

                dids = [tup[0] for tup in self.sql.fetchall()]

                for did in dids:

                    # for each drug, search for a connection

                    results = PGKB_connect(authobj, var, did)

                    if results is not None:

                        for phen in results:

                            if var_alleles[var] in phen['allele']:

                                self.advice.setdefault(did, []).append((var, phen['phenotype']))

tom = Patient('data/hpc/STE0097_WGS_PharmacoGenomics.vcf')
app = Interpreter(tom)
app.DrugAdvice()
print app.advice