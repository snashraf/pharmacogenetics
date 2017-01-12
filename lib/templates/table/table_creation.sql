CREATE TABLE Pairs (
	PairID binary PRIMARY KEY AUTOINCREMENT,
	GeneID text,
	DrugID text,
	GuideID text
);

CREATE TABLE Drugs (
	DrugID text,
	ChemName text
);

CREATE TABLE Genes (
	GeneID text,
	GeneSymbol text,
	Chromosome text,
	Start integer,
	End integer
);

CREATE TABLE DrugNames (
	DrugID text,
	TradeName text
);

CREATE TABLE DrugTerms (
	DrugID text,
	DrugTerm text
);

CREATE TABLE Haplotypes (
	HapID text,
	GeneID text,
	Starname text,
	HGVS text
);

CREATE TABLE HapVars (
	HapID text,
	RSID text,
	AltAllele text,
	MutType text
);

CREATE TABLE Variants (
	VarID text,
	RSID text,
	Description text
);

CREATE TABLE LocPGKB (
	VarID text,
	Chromosome text,
	Start integer,
	End integer,
	RefAllele text
);

CREATE TABLE LocVCF (
	VarID text,
	Chromosome text,
	Start integer,
	End integer,
	RefAllele text
);

CREATE TABLE AltAlleles (
	VarID text,
	AltPGKB text,
	AltVCF text
);

CREATE TABLE AnnoVar (
	AnID text,
	PairID binary,
	VarID text,
	LevelOfEvidence text
);

CREATE TABLE AnnoHap (
	AnID text,
	PairID binary,
	HapID text,
	LevelOfEvidence text
);

CREATE TABLE AnnoText (
	AnID text,
	Allele text,
	Phenotype text
);

CREATE TABLE GuideOptions (
	GuideID text,
	HapID binary,
	HapName text
);

