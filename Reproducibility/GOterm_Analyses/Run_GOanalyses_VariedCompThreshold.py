
from goatools import obo_parser
from goatools.go_enrichment import GOEnrichmentStudy
from Bio import SeqIO

def main():

    go_obo = 'go-basic.obo'
    go = obo_parser.GODag(go_obo)

    proteome = 'UP000005640_9606'
    abbrev = 'Hsapiens'
    
    filenames = ['TableSX_NewSR-related_RBPs_60percCompThreshold_GOresults', 'TableSY_NewSR-related_RBPs_65percCompThreshold_GOresults']
    prot_df = get_hits()
    
    assoc = get_assoc(abbrev)
    i = 0
    for comp_threshold in ['60', '65']:
        filename = filenames[i]
        lcd_prots = prot_df[comp_threshold]
        background_prots = get_background_prots(proteome)
        
        if len(lcd_prots) > 0:
            methods = ["bonferroni", "sidak", "holm", "fdr"]
            g = GOEnrichmentStudy(background_prots, assoc, go, propagate_counts=True, alpha=0.05)
            g_res = g.run_study(lcd_prots)
            filename += '.tsv' #change file type here to txt or xlsx if desired
            g.wr_tsv(filename, g_res)
            
        i += 1


def get_hits():
    h = open('Table S5 - Human SR and SR-related Proteins by S-R Composition Threshold.tsv')
    header = h.readline()
    
    df = {}
    for line in h:
        items = line.rstrip().split('\t')
        comp_threshold = items[0]
        sr_and_rbp_prots = items[8].split(', ')
        df[comp_threshold] = sr_and_rbp_prots
    h.close()
    
    return df

        
def get_background_prots(proteome):
    
    prots = []
    h = open(proteome + '.fasta')
    for seq_record in SeqIO.parse(h, 'fasta'):
        id = str(seq_record.id)
        junk, uniprot_id, junk = id.split('|')
        prots.append( uniprot_id )
    h.close()
        
    return prots
    
    
def get_assoc(abbrev):
    """
    Reads the proper gene association file (GAF).
    Return is a dictionary with Protein IDs as keys, and a set of associated GO terms as values (this is a formal Python set, not a list)
    
    Each list represents a separate GO term associated with the gene key, and contains the following:
    0) Database ("UniProtKB")
    1) Protein ID
    2) Common gene name
    3) Qualifier (e.g. 'NOT', 'contributes_to', 'colocalizes_with', etc. This is optional and is sometimes an empty string)
    4) GO_ID (i.e. the GO term associated with the gene. NOTE: multiple GO terms associated with a single gene will be in separate entries in the dictionary)
    5) Literature reference that the annotation was derived from
    6) Evidence code
    7) Mystery - don't really know what this column represents and doesn't seem to match README file description on GO website
    8) Aspect (i.e. molecular function, biological process, or cellular component...indicated with abbreviation)
    9) Object symbol (e.g. 'Mitochondrial 21s RNA)
    10) Object synonym
    11) Object type (i.e. gene, protein, etc.)
    12) Taxon ID (559292 for s cerevisiae)
    13) Date of annotation
    14) Assigned by (a database)
    """

    h = open(abbrev + '.gaf')
    
    all_gos = []
    assoc = {}
    
    #creates 'pop' list and 'assoc' dictionary to pass to GOEnrichmentStudy() module in goatools
    for line in h:
        if line.startswith('!'):
            continue
        items = line.rstrip().split('\t')
        temps = items[10].split('|')
        uniprot_id = items[1]
        assoc[uniprot_id] = assoc.get(uniprot_id, [])
        assoc[uniprot_id].append(items[4])
        all_gos.append(items[4])
        
    set_gos = set(all_gos)

    for key in assoc:
        assoc[key] = set(assoc[key])
    
    h.close()

    return assoc


if __name__ == '__main__':
    main()