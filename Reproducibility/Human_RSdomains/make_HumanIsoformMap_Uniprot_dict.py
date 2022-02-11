
import pickle
from Bio import SeqIO

def main():

    df, iso_to_orig = get_orig_prots()
    h = open('UP000005640_9606_additional.fasta')
    
    for seq_record in SeqIO.parse(h, 'fasta'):
        desc = str(seq_record.description)
        junk, uniprot, *junk = desc.split('|')
        junk, temp = desc.split('Isoform of ')
        isoform_of, *junk = temp.split(', ')
        
        df[isoform_of].append( uniprot )
        iso_to_orig[uniprot] = isoform_of
        
    h.close()
    
    pf = open('Uniprot_OrigProt_to_IsoformList_MappingDict.dat', 'wb')
    pickle.dump(df, pf)
    pf.close()
    
    pf = open('Uniprot_Isoform_to_OrigProt_OneToOne_MappingDict.dat', 'wb')
    pickle.dump(iso_to_orig, pf)
    pf.close()
    
    
def get_orig_prots():

    h = open('UP000005640_9606.fasta')
    df = {}
    iso_to_orig = {}
    for seq_record in SeqIO.parse(h, 'fasta'):
        id = str(seq_record.id)
        junk, uniprot, *junk = id.split('|')
        df[uniprot] = []
        iso_to_orig[uniprot] = uniprot
        
    h.close()
        
    return df, iso_to_orig


if __name__ == '__main__':
    main()