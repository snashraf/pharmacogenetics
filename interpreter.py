#!/usr/bin/python
# -*- coding: utf-8 -*-
from patient import Patient
from data import DataCollector
import json
import sqlite3
import pprint
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
    try:
        r = client.get(uri)
    except:
        return None

    data = r.json()

    # if there is just a dict and not a list it is an error message

    if type(data) == dict:
        if "errors" not in str(data):
            return data
        else:
            return None

    # otherwise return valid results

    elif type(data) == list:
        return data


def PGKB_connect(authobj, mode, a, did):
    """
    Find connections between two pharmgkb objects.
    Usually a haplotype or rsid compared with a drug id.
    :param authobj: authentication session resulting from Authenticate()
    :param a: object 1 to compare
    :param did: object 2 to compare, usually a drug id
    :return:
    """
    if mode == "clinAnno":
        uri = \
            'https://api.pharmgkb.org/v1/report/pair/%s/%s/clinicalAnnotation' \
            % (a, did)

    # get json file from filled in uri, use authentication client

    result = getJson(uri, authobj)
    results = []

    # collect matching results
    if result is not None:
        for doc in result:
            if mode == "clinAnno":

                for phen in doc['allelePhenotypes']:
                    results.append(phen)

            elif mode =="clinGuide":

                for item in doc:
                    print doc
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
        #self.DrugAnnotations()
        self.conn.commit()

    def DrugAnnotations(self):

        self.sql.execute("DROP TABLE IF EXISTS patientadvice")
        self.sql.execute('''CREATE TABLE patientadvice
                                            (did text, gid text, varid text, allele text, advice text, dose text)''')
        authobj = Authenticate()

        self.hapdrug = {}

        # from haplotypes = [PA..., PA...]
        # from rsids = [rs..., rs..., rs...]

        self.sql.execute('''SELECT d.rsid, p.chr1, p.chr2
                            FROM patientvars p
                            JOIN design d ON d.pos = p.pos
                            JOIN variants v ON d.rsid = v.rsid''')

        var_alleles = {rsid:( chr1 + chr2 )
                    for( rsid, chr1, chr2 )
                    in self.sql.fetchall() }

        self.sql.execute("SELECT DISTINCT hapid FROM patienthaps")

        self.haplotypes = [tup[0] for tup in self.sql.fetchall()]

        input = self.haplotypes + var_alleles.keys()

        for var in input:

            # check for haplotype or rs nomenclature

            if "PA" in var:
                self.sql.execute("SELECT starhap FROM patienthaps WHERE hapid = ?", (var,))
                var_alleles[var] = self.sql.fetchone()[0]
                self.sql.execute("SELECT DISTINCT gid FROM alleles WHERE hapid = ?", (var,))

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

                    results = PGKB_connect(authobj, "clinAnno", var, did)

                    if results is not None:

                        for phen in results:

                            allele = var_alleles[var]

                            if allele in phen['allele']:

                                advice = phen['phenotype']

                                item = (did, gid, var, allele, advice, "N/A") # tuple with 5 items

                                self.sql.execute('''INSERT INTO patientadvice VALUES(?,?,?,?,?,?)''', item)

    def DoseGuide(self):
        authobj = Authenticate()
        # input: haplotypes and drug matches
        self.sql.execute('SELECT gid,did,allele FROM patientadvice WHERE varid LIKE "%PA%"')
        results = self.sql.fetchall()
        sources = ["cpic", "dpwg", "pro"]
        for (gid, did, allele) in results:
            for source in sources:
                uri = "https://api.pharmgkb.org/v1/data/guideline?source=%s&relatedChemicals.accessionId=%s&relatedGenes.accessionId=%s&view=max" \
                        % (source, did, gid)
                guidelines = getJson(uri, authobj)
                if guidelines is not None:
                    for doc in guidelines:
                        if doc["objCls"] == "Guideline":
                            guid = doc["id"]
                        else:
                            break
                    uri = "https://api.pharmgkb.org/v1/report/guideline/%s/options" \
                        % guid
                    doseOptions = getJson(uri, authobj)
                    uri = 'https://api.pharmgkb.org/v1/data/guideline/%s?view=max'\
                        % guid
                    doseInfo = getJson(uri, authobj)
                    pp = pprint.PrettyPrinter()
                    pp.pprint(doseInfo)
                    if len(doseOptions["data"]) > 0:
                        for choice in doseOptions["data"]:
                            gene = choice['symbol']
                            options = choice['options']
                            score = 0
                            haps = allele.split("/")
                            for chr in haps:
                                if chr in options:
                                    score += 1
                            if score == 2:
                                uri = "https://api.pharmgkb.org/v1/report/guideline/%s/annotations?genes=%s&alleles=%s&alleles=%s" \
                                    % (guid, gene, haps[0], haps[1])
                                doseGuide = getJson(uri, authobj)
                                if doseGuide is not None:
                                    print doseGuide

tom = Patient('data/hpc/test.vcf')
app = Interpreter(tom)
app.DoseGuide()

"""
NOTE TO SELF: USE THIS

steps to take:
HAPLOTYPE PROPERLY ;W;
find the remaining variants... check the github repos i dled for some tables
and convert the rest back, hopefully

find a guideline for gene-drug pairsA
https://api.pharmgkb.org/v1/data/guideline?source=cpic&relatedChemicals.accessionId=PA449088&relatedGenes.accessionId=PA128&view=max

get options:
https://api.pharmgkb.org/v1/report/guideline/PA166114461/options

find match and corresponding advice in options (after haplotype assignment)
https://api.pharmgkb.org/v1/report/guideline/PA166114461/annotations?genes=CFTR&alleles=G1244E&alleles=G1244E

tadaa!

also check these people did something similar:
http://www.sustc-genome.org.cn/vp/
"""