
import pickle

def main():

    uniprot_PTMs = {}
    h = open('ActiveDriverDB2021_PTMsites.tsv')
    header = h.readline()
    for line in h:
        items = line.rstrip().split('\t')
        position = int(items[1])
        res = items[2]
        types = items[3].split(',')
        uniprot = items[7]

        for mod_type in types:
            if 'SARS-CoV-2' in mod_type:    #FILTER OUT SITES THAT ARE SPECIFICALLY OBSERVED DURING SARS-COV-2 INFECTION: THESE MAY BE NON-NATIVE SITES. MOST OF THESE SITES ARE STILL INCLUDED IN THE DATASET BECAUSE THE PTM IS ALSO OBSERVED IN THE ABSENCE OF SARS-COV-2 INFECTION.
                continue
            uniprot_PTMs[uniprot] = uniprot_PTMs.get(uniprot, [])
            uniprot_PTMs[uniprot].append( (res, position, mod_type) )
            
    h.close()

    pf = open('ActiveDriverDB2021_PTMs_mapped_to_UniprotIDs.dat', 'wb')
    pickle.dump(uniprot_PTMs, pf)
    pf.close()


if __name__ == '__main__':
    main()