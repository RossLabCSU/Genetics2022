
import pickle
from Bio import SeqIO
pfam_df = pickle.load(open('Human_Pfam_Annotations_by_UniprotAccession.dat', 'rb'))
gene_conversion = pickle.load(open('EnsemblGeneID_to_UniProtID_converter.dat', 'rb'))

def main():

    pfam_rrm_accs = get_RRM_pfam_accs()
    gerstberger_rbps = get_gerstberger_rbps()
    goterm_rbps = get_goterm_rbps()
    
    output = open('Complete_list_of_RBPs.txt', 'w')
    
    track_prots = set()
    
    h = open('UP000005640_9606.fasta')
    for seq_record in SeqIO.parse(h, 'fasta'):
        id = str(seq_record.id)
        junk, uniprot, *junk = id.split('|')
        
        # CHECK FOR RRM
        contains_rrm = 0
        if uniprot in pfam_df:
            for annot in pfam_df[uniprot]:
                pfam_acc = annot[4]
                if pfam_acc in pfam_rrm_accs:
                    contains_rrm = 1
                    
        # CHECK FOR OVERLAP WITH GERSTBERGER DATASET
        is_gerstberger_rbp = 0
        if uniprot in gerstberger_rbps:
            is_gerstberger_rbp = 1
            
        # CHECK FOR ANNOTATION WITH "RNA Binding" GO TERM
        is_goterm_rbp = 0
        if uniprot in goterm_rbps:
            is_goterm_rbp = 1
            
        if sum([contains_rrm, is_gerstberger_rbp, is_goterm_rbp]) > 0 and uniprot not in track_prots:
            output.write(uniprot + '\n')
            track_prots.add(uniprot)
            
    output.close()
            
    
def get_RRM_pfam_accs():

    h = open('Pfam_Accessions_for_RRMclan_Only.tsv')
    pfam_rrm_accs = []
    for line in h:
        items = line.rstrip().split('\t')
        pfam_acc = items[0]
        pfam_rrm_accs.append(pfam_acc)
        
    h.close()
    
    return pfam_rrm_accs
    
                
def get_gerstberger_rbps():

    h = open('RBPs_TableS1_from_Gerstberger2014.txt')
    header = h.readline()

    rbps = []
    for line in h:
        items = line.rstrip().split('\t')
        ensembl = items[2]
        
        if ensembl in gene_conversion:
            for prot in gene_conversion[ensembl]:
                rbps.append(prot)
            
    h.close()

    return rbps
    

def get_goterm_rbps():

    h = open('Hsapiens_RNAbinding_Proteins_GOid0003723.txt')
    prots = []
    for line in h:
        prots.append(line.rstrip())
    h.close()
        
    return prots
        

if __name__ == '__main__':
    main()