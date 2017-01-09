#!/usr/bin/python
# -*- coding: utf-8 -*-

import sqlite3

# --------------------------------------------------------------------------

class Database(object):

    '''
....Database class. Initializes database and takes care of doing big changes.

....'''

    def __init__(self, dbname):

        self.conn = sqlite3.connect('%s.db' % dbname)  # connect to db- if it doesn't exist, create it

        self.sql = self.conn.cursor()  # cursor for sqlite3, used to do things in database

        self.setDefaults()

    def setDefaults(self):

        self.tableoptions = {

            'drugpairs': 'gid text, did text, starhaps text, options text, rsids text',
        
            'genes': 'gid text UNIQUE, symbol text, chr text, start text, stop text',
        
            'alleles': "hapid text, gid text, starname text,\
						hgvs text, rsid text, alt text,\
						UNIQUE(hapid, rsid, starname)\
						ON CONFLICT REPLACE)",
        
            'variants': "rsid text, varid text, gid text, chr text,\
					pbegin int, pend int, pref text, palt text,\
					vbegin int, vend int, vref text, valt text,\
					type text,\
									UNIQUE(rsid,gid)\
									ON CONFLICT REPLACE)",
        
            'alias': 'rsid text, alias varchar(255) PRIMARY KEY',
        
            'drugs': 'did text, name text, terms text',
        
            }


    def removeTable(self, tabname):
        """
........This function rebuilds the database, use for updating the db periodically?
........:return:
........"""

        self.sql.execute("DROP TABLE IF EXISTS {}".format(tabname))


    def alterDefault(self, tabname, tabvalues):

        self.tableoptions[tabname] = tabvalues


    def createTable(self, tabname):

        self.sql.execute('CREATE TABLE {}({})'.format(tabname, self.tableoptions[tabname]))
        
        self.conn.commit()


    def insertValues(self, tabname, tabvalues):

    	qmarks = ",".join(["?" for item in tabvalues])

        self.sql.execute('INSERT INTO {} VALUES({})'.format(tabname, qmarks), (tabvalues))

        self.conn.commit()


    def remakeTable(self, tabname):

    	self.removeTable(tabname)

    	self.createTable(tabname)

    	self.conn.commit()