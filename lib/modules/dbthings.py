#!/usr/bin/python
# -*- coding: utf-8 -*-

import sqlite3

conn = sqlite3.connect('pharmacogenetics.db')
sql = conn.cursor()

sql.execute('SELECT DISTINCT rsid, alt from alleles;')
mutations = sql.fetchall()

headers = {
    'CYP1A2': 'NM_000761.3',
    'CYP2A6': 'NM_000762.5',
    'CYP2B6': 'NM_000767.4',
    'CYP2C8': 'NM_000770.3',
    'CYP2C9': 'NG_008385.1',
    'CYP2D6': 'M33388',
    'CYP3A4': 'NM_001202855.2',
    'UGT1A1': 'AF297093.1',
    'NAT2': 'X14672',
    }

middle = ':c.'

sql.execute('DROP TABLE IF EXISTS newvariants')
sql.execute('CREATE TABLE newvariants(name text, header text, pos text, ref text, alt text)'
            )

prev_val = ''

for (rsid, alt) in mutations:
    sql.execute("select alt from alleles where hgvs like '%=%' and rsid = ?"
                , (rsid, ))

    try:
        ref = sql.fetchone()[0]
    except:
        ref = 'NA'

    if 'rs' in rsid and len(rsid) < 15:
        if ref == alt:
            continue
        elif ref != alt:
            item = (rsid, 'n/a', 'n/a', ref, alt)
    elif 'rs' not in rsid:

        if ':' in rsid and '(' not in rsid:
            symbol = rsid.split(':')[0]
            if symbol not in headers.keys():
                seq = symbol
            else:
                seq = headers[symbol]
        else:
            symbol = rsid.split(' ')[0]
            seq = headers[symbol]

        pos = ''.join(rsid.split(' ')[1::])
        item = (rsid, seq, pos, ref, alt)

    sql.execute('INSERT INTO newvariants VALUES(?,?,?,?,?)', item)

conn.commit()