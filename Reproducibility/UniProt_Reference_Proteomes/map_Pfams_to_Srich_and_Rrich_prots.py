
import pickle
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']

def main():

    residues = ['R', 'S']
    for i in range(len(residues)):
        aa = residues[i]
        if aa == 'R':
            limited_aa = 'S'
        else:
            limited_aa = 'R'
            
        for domain in domains:
            annots = pickle.load(open(domain + '_' + aa + 'richOnly_Gathered_Pfam_Results.dat', 'rb'))
            filename = domain + '_' + aa + '-rich-only_RESULTS'
            h = open(filename + '.tsv')
            for i in range(7):
                h.readline()
            output = open(domain + '_' + aa + '-richOnly_WithPfamAnnotations.tsv', 'w')
            
            header = h.readline().rstrip().split('\t')
            output.write('\t'.join( header[:7] + ['Pfam Domain IDs', 'Pfam Domain Names', 'Pfam Domain Boundaries', 'Pfam E-values']) + '\n')

            pfam_ids = set()
            id_list = []
            
            for line in h:
                items = line.rstrip().split('\t')
                uniprot = items[2]
                junk, uniprot, *junk = uniprot.split('|')
                if uniprot in annots:
                    pfams = annots[uniprot]
                    accs = []
                    names = []
                    descs = []
                    bounds = []
                    bitscores = []
                    evals = []
                    for pfam in pfams:
                        env_start = pfam['env']['from']
                        env_end = pfam['env']['to']
                        pfam_name = pfam['name']
                        pfam_acc = pfam['acc']
                        pfam_desc = pfam['desc']
                        pfam_eval = pfam['evalue']
                        pfam_bitscore = pfam['bits']
                        pfam_clan = pfam['clan']
                        
                        accs.append(pfam_acc)
                        descs.append(pfam_desc)
                        bounds.append( '(' + env_start + '-' + env_end + ')' )
                        names.append(pfam_name)
                        bitscores.append(bitscores)
                        evals.append(pfam_eval)

                    output.write('\t'.join( items[:7] + [';'.join(accs), ';'.join(names), ';'.join(bounds), ';'.join(evals)] ) + '\n')

            h.close()
            output.close()


if __name__ == '__main__':
    main()