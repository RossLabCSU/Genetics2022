
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
            if 'SARS-CoV-2' in mod_type:
                continue
            uniprot_PTMs[uniprot] = uniprot_PTMs.get(uniprot, [])
            uniprot_PTMs[uniprot].append( (res, position, mod_type) )
            
    h.close()

    pf = open('ActiveDriverDB2021_PTMs_mapped_to_UniprotIDs.dat', 'wb')
    pickle.dump(uniprot_PTMs, pf)
    pf.close()


if __name__ == '__main__':
    main()