#!/usr/bin/python
# -*- coding: utf-8 -*-

import sqlite3
import os
import re
import glob
from jinja2 import Template, FileSystemLoader, Environment

from modules.pgkb_functions import *

# --------------------------------------------------------------------------

class Database(object):

    '''
 Database class. Initializes database and takes care of doing big changes.

 '''

    def __init__(self, dbname):

        self.conn = sqlite3.connect('%s.db' % dbname)  # connect to db- if it doesn't exist, create it

        self.sql = self.conn.cursor()  # cursor for sqlite3, used to do things in database

        self.path = os.getcwd()

        self.tempfolder = self.path + '\\lib\\templates'

        self.setDefaults()

    def loadSQL(self, path):

        conv = ''
        filt = '\t\n'

        with open(path, 'r') as f:
            for line in f.readlines():
                conv += re.sub('[\n\t]', '', line)

        return conv


    def insertSQL(self, tabname):

		templateLoader = FileSystemLoader( searchpath=self.tempfolder)
		
		templateEnv = Environment( loader=templateLoader )
		
		TEMPLATE_FILE = tabname + '.ins'
		
		template = templateEnv.get_template( TEMPLATE_FILE )

		return template


    def setDefaults(self):

        dropTables = glob.glob(self.tempfolder + "*.rm")
        createTables = glob.glob(self.tempfolder + "*.tab")

        for template in dropTables + createTables:

            sql = self.loadSQL(template)

            self.sql.executescript(sql)

        self.conn.commit()


    def removeTable(self, tabname):

        sql = self.loadSQL(self.tempfolder + '\\' + tabname + '.rm')

        self.sql.executescript(sql)

        self.conn.commit()


    def createTable(self, tabname):

        sql = self.loadSQL(self.tempfolder + '\\' + tabname + '.tab')

        self.sql.executescript(sql)

        self.conn.commit()

    def remakeTable(self, tabname):

        self.removeTable(tabname)

        self.createTable(tabname)

        self.conn.commit()

# -----------------------------------------------------

db = Database('test')