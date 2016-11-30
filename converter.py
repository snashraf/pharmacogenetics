import sqlite3

conn = sqlite3.connect('pharmacogenetics.db')
sql = conn.cursor()

sql.execute("SELECT DISTINCT rsid, alt from alleles ;")
