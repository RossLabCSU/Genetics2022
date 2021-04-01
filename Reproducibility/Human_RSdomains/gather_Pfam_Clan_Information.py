
import pickle

def main():

    h = open('Pfam-A.clans.tsv')
    
    df = {'DESCRIPTION':'\nDictionary containing Pfam domain accessions as keys, and value is a list containing the associated Pfam domain information.\nItems in the list are in the following index order:\n0)Pfam domain accession\n1)Pfam clan associated with Pfam domain accession\n2)Clan short name?\n3)Pfam short name?\n4)Pfam domain description\n'}
    output = open('Pfam_Accessions_for_RRMclan_Only.tsv', 'w')

    rrm = []
    for line in h:
        items = line.rstrip().split('\t')
        pfam_acc = items[0]
        clan_info = items[2]
        df[pfam_acc] = items[1:]
        if clan_info == 'RRM':
            output.write(line)
            rrm.append(items)

    output.close()
    h.close()

    pf = open('Pfam_Info_by_PfamAccession.dat', 'wb')
    pickle.dump(df, pf)
    pf.close()

if __name__ == '__main__':
    main()