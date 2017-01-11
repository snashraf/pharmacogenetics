#!/usr/bin/python
# -*- coding: utf-8 -*-

import sqlite3
from jinja2 import Template
import re
import os
# --------------------------------------------------------------------------

class Database(object):

    '''
	Database class. Initializes database and takes care of doing big changes.

	'''

    def __init__(self, dbname):

        self.conn = sqlite3.connect('%s.db' % dbname)  # connect to db- if it doesn't exist, create it

        self.sql = self.conn.cursor()  # cursor for sqlite3, used to do things in database

        self.path = os.getcwd()

        self.tempfolder = self.path + "\\templates\\table"

        self.insfolder = self.path + "\\templates\\insert"

        self.setDefaults()

    def formatSQL(self, path):
    
        conv = ""
        filt = '\t\n'
        
        with open(path, "r") as f:
            for line in f.readlines():
                conv += re.sub('[\n\t]', '', line)

        return conv
        

    def renderSQL(self, sql_string, context):

        template = Template(sql_string)

        sql = template.render(json = context)

        return sql


    def insertSQL(self, tabname, context):

        conv = self.formatSQL(self.insfolder + "\\" + tabname + ".txt")

        sql = self.renderSQL(conv, context)

        self.sql.execute(sql)


    def setDefaults(self):

        templates = os.listdir(self.tempfolder)

        for template in templates:

            self.removeTable(template.rstrip(".txt"))

            sql = self.formatSQL(self.tempfolder + "\\" + template)

            self.sql.execute(sql)

        self.conn.commit()


    def removeTable(self, tabname):
        """
		This function rebuilds the database, use for updating the db periodically?
		return:
		"""

        self.sql.execute("DROP TABLE IF EXISTS {}".format(tabname))


    def alterDefault(self, tabname, tabvalues):

        self.tableoptions[tabname] = tabvalues


    def createTable(self, tabname):

    	sql = self.formatSQL(self.tempfolder + "\\" + tabname + ".txt")

        self.sql.execute(sql)

        self.conn.commit()


    def insertValues(self, tabname, tabvalues):

    	qmarks = ",".join(["?" for item in tabvalues])

        self.sql.execute('INSERT INTO {} VALUES({})'.format(tabname, qmarks), (tabvalues))


    def remakeTable(self, tabname):

    	self.removeTable(tabname)

    	self.createTable(tabname)

    	self.conn.commit()

db = Database("test")