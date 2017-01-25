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

		print path

		self.tempfolder = path.replace("db", "templates")

		self.conn = sqlite3.connect(dbpath)  # connect to db- if it doesn't exist, create it

		self.sql = self.conn.cursor()  # cursor for sqlite3, used to do things in database

		templateLoader = FileSystemLoader( searchpath=self.tempfolder)

		self.templateEnv = Environment( loader=templateLoader )

	def loadSQL(self, path):

		sql = ''

		with open(path, 'r') as f:

			for line in f.readlines():

				sql += re.sub('[\n\t]', '', line)

		return sql

	def getTemplate(self, template):

		TEMPLATE_FILE = template

		template = self.templateEnv.get_template( TEMPLATE_FILE )

		return template

	def insertSQL(self, tabname):

		TEMPLATE_FILE = tabname + '.ins'

		template = self.templateEnv.get_template( TEMPLATE_FILE )

		return template

	def setDefaults(self):

		dropTables = glob.glob(os.path.join(self.tempfolder, "*.rm"))

		createTables = glob.glob(os.path.join(self.tempfolder, "*.tab"))

		for template in dropTables + createTables:

			sql = self.loadSQL(template)

		self.sql.executescript(sql)

		self.conn.commit()

		# test

	def removeTable(self, tabnames):
		if type(tabnames) != list:
			tabnames = []
		for tabname in tabnames:
			sql = self.loadSQL(os.path.join(self.tempfolder, tabname + '.rm'))
			self.sql.executescript(sql)
		self.conn.commit()

	def createTable(self, tabnames):
		if type(tabnames) != list:
			tabnames = []
		for tabname in tabnames:
			sql = self.loadSQL(os.path.join(self.tempfolder, tabname + '.tab'))
			self.sql.executescript(sql)
		self.conn.commit()


	def remakeTable(self, tabname):

		self.removeTable(tabname)

		self.createTable(tabname)

		self.conn.commit()

# -----------------------------------------------------
