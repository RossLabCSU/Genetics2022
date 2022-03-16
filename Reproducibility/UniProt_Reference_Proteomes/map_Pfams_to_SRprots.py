
import pickle
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
file_leaders = {'Archaea':'TableS8', 'Bacteria':'TableS9', 'Eukaryota':'TableS10', 'Viruses':'TableS11'}

def main():

    for domain in domains:
        annots = pickle.load(open(domain + '_Gathered_Pfam_Results.dat', 'rb'))
        filename = domain + '_SR_proteins_with_Combined_S-R_Above_70_MERGED_DOMAINS'
        h = open(filename + '.tsv')
        output = open(file_leaders[domain] + '_' + filename + '_WithPfamAnnotations.tsv', 'w')
        
        header = h.readline().rstrip().split('\t')
        output.write('\t'.join( header + ['Pfam Domain IDs', 'Pfam Domain Names', 'Pfam Domain Boundaries', 'Pfam E-values']) + '\n')

        pfam_ids = set()
        id_list = []
        
        for line in h:
            items = line.rstrip().split('\t')
            uniprot = items[2]
            accs = ['N/A']
            names = ['N/A']
            descs = ['N/A']
            bounds = ['N/A']
            bitscores = ['N/A']
            evals = ['N/A']
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

            output.write('\t'.join( items + [';'.join(accs), ';'.join(names), ';'.join(bounds), ';'.join(evals)] ) + '\n')

        h.close()
        output.close()


if __name__ == '__main__':
    main()
