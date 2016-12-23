from optparse import OptionParser
from data import DataCollector

def main():
    
    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")
    
    parser.add_option("-d", "--makedb",
                      action="append",
                      dest="tables",
                      default=[],
                      help="Download database (specify tables, default is do whole database)\nOptions are pairs (gene-drug pairs), genes, vars (variants) and drugs(chemicals)")
    
    parser.add_option("-p", "--patient",
                      action="store", # optional because action defaults to "store"
                      dest="gvcf",
                      default=None,
                      help="Patient g.vcf file to parse",)
    
    (options, args) = parser.parse_args()
        
    if len(options.tables) > 0:
        
        CreateDB(options.tables)
    
    if options.gvcf:
        
        if ".gz" not in options.gvcf:
            
            print "Please convert to .gz and create tabix file first."
        
        else:
            
            pass
    
    
def CreateDB(tables):
        
    d = DataCollector()

    if "all" in tables:
                    
        d.Update()
        
    else:
        
        for table in tables:
            
            if table == "pairs":
                
                d.GetPairs()
            
            elif table == "genes":
              
                d.GetGeneData()
            
            elif table == "allvars":
            
                d.GetVarData()

		d.GetNonRS()

	    elif table == "rsvars":

		d.GetVarData()
            
            elif table == "drugs":
                
                d.GetChemData()

	    elif table == "hapvars":

		d.GetNonRS()

        d.conn.commit()

if __name__ == '__main__':
    main()
