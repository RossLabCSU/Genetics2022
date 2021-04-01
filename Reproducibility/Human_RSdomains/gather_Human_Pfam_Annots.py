
import pickle

def main():

    h = open('9606_Pfam_Domains_Mapping_Uniprot.tsv')
    for i in range(3):
        h.readline()
    
    # If trying to access the Pfam short domain name, use index 5 for each sublist
    df = {'DESCRIPTION':'\nDictionary containing Uniprot accessions for human proteins as keys, and values are a list of lists where each sublist contains the Pfam domain information.\nItems in the list are in the following index order:\n0)Pfam domain start site\n1)Pfam domain end site\n2)Envelope start\n3)Envelope end\n4)Pfam accession\n5)Pfam short name\n6)Domain type\n7)HMM start\n8)HMM end\n9)HMM length\n10)Bit score\n11)E-value\n12)Clan\n'}
    for line in h:
        items = line.rstrip().split('\t')
        prot = items[0]
        info = items[1:]
        df[prot] = df.get(prot, [])
        df[prot].append(info)
        
    h.close()

    pf = open('Human_Pfam_Annotations_by_UniprotAccession.dat', 'wb')
    pickle.dump(df, pf)
    pf.close()  

if __name__ == '__main__':
    main()