#!/usr/bin/python
# -*- coding: utf-8 -*-

import sqlite3


class DocHandler:

    def __init__(self):

        self.conn = sqlite3.connect('pharmacogenetics.db')

        self.sql = self.conn.cursor()

    def AnnOverview(self):

            # methods: basic text document, latex

            # group by: drug | gene

        self.sql.execute('select distinct a.did, c.name from annotations a join chemicals c on a.did = c.did'
                         )

        hline = 35 * '-'

        print 'Overview Drug Annotations'

        print hline

        lst = self.sql.fetchall()

        for (did, name) in lst:

            self.sql.execute('select terms from chemicals where did = ?'
                             , (did, ))

            terms = [tup[0] for tup in self.sql.fetchall()]

            print hline

            print name.capitalize(), '(%s)' % ', '.join(terms)

            print hline

            self.sql.execute('select a.advice, g.symbol, a.loe from annotations a join variants v on a.varid=v.varid join genes g on v.gid=g.gid where a.did = ? order by loe asc'
                             , (did, ))

            advices = self.sql.fetchall()

            for (advice, symbol, loe) in advices:

                print '|| GENE: %s LevelOfEvidence: %s ||' % (symbol,
                        loe), advice, '\n'


d = DocHandler()
d.AnnOverview()