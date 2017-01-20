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

    def __init__(self, dbpath):

        path = os.path.dirname(dbpath)

        self.tempfolder = path.replace("db", "templates")

        # connect to db- if it doesn't exist, create it
        self.conn = sqlite3.connect(dbpath)

        self.sql = self.conn.cursor()  # cursor for sqlite3, used to do things in database

    def loadSQL(self, path):

        sql = ''

        with open(path, 'r') as f:

            for line in f.readlines():

                sql += re.sub('[\n\t]', '', line)

        return sql

    def insertSQL(self, tabname):

        templateLoader = FileSystemLoader(searchpath=self.tempfolder)

        templateEnv = Environment(loader=templateLoader)

        TEMPLATE_FILE = tabname + '.ins'

        template = templateEnv.get_template(TEMPLATE_FILE)

        return template

    def setDefaults(self):

        dropTables = glob.glob(os.path.join(self.tempfolder, "*.rm"))

        createTables = glob.glob(os.path.join(self.tempfolder, "*.tab"))

        for template in dropTables + createTables:

            sql = self.loadSQL(template)

        self.sql.executescript(sql)

        self.conn.commit()

        # test

    def removeTable(self, tabname):

        sql = self.loadSQL(os.path.join(self.tempfolder, tabname + '.rm'))

        self.sql.executescript(sql)

        self.conn.commit()

    def createTable(self, tabname):

        sql = self.loadSQL(os.path.join(self.tempfolder, tabname + '.tab'))

        self.sql.executescript(sql)

        self.conn.commit()

    def remakeTable(self, tabname):

        self.removeTable(tabname)

        self.createTable(tabname)

        self.conn.commit()

# -----------------------------------------------------
