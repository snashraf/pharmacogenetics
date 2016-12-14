import sqlite3



conn = sqlite3.connect('pharmacogenetics.db')

sql = conn.cursor()

sql.execute("DROP TABLE IF EXISTS patienthaps")

sql.execute("CREATE TABLE IF NOT EXISTS patienthaps(gid text, hapid text, score int)")

sql.execute("SELECT DISTINCT gid from genes")

gids = [tup[0] for tup in sql.fetchall()]

# get list of all gids

for gid in gids:

	sql.execute("SELECT rsid, alt from alleles where hgvs like '%=%' and rsid like '%rs%' and gid=?", (gid,))

	reference = { rsid : alt for (rsid, alt) in sql.fetchall() }
	
	# get list of all hapids for this gene
	
	sql.execute("SELECT DISTINCT hapid, starname, hgvs from alleles where gid=?", (gid,))
	
	hapids = sql.fetchall()
	
	for (hapid, starhap, hgvs) in hapids:
		
		if "[0]" not in hgvs and "=" not in hgvs:
			
			# get haplotype alleles and create complete dictionary
			
			sql.execute("SELECT rsid, alt from alleles where hapid=? and rsid like '%rs%'", (hapid,))

			haprsids = { rsid : alt for (rsid, alt) in sql.fetchall()}

			haprsids = dict(reference, **haprsids)
			
			# get patient rsids
			
			sql.execute("select distinct v.rsid, p.alt from variants v join patientvars p on p.start=v.start join alleles a on v.gid=a.gid where p.alt != '.' and a.hapid = ?", (hapid,))
			
			patrsids = { rsid : alt for (rsid, alt) in sql.fetchall()}
			
			patrsids = dict(reference, **patrsids)

			# compare patient to haplotype rsids
			
			modified = {k : (patrsids[k], haprsids[k]) for k in patrsids.keys() if patrsids[k] != haprsids[k]}
			
			score = len(modified)

# -------------------------------------------------------------------------------------------------------------------------------------------------------

			item = (gid, hapid, score)
			
			sql.execute("INSERT INTO patienthaps VALUES(?,?,?)", item)

# -------------------------------------------------------------------------------------------------------------------------------------------------------

conn.commit()
