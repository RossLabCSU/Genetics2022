
file_leaders = {'Archaea':'TableS6', 'Bacteria':'TableS7', 'Eukaryota':'TableS8', 'Viruses':'TableS9'}

def main(args):

    threshold = 70
    domain = args.domain

    df = {}
    sr_prots = {}
    for s_comp in range(20, 105, 5):
        for r_comp in range(20, 105, 5):
            if s_comp + r_comp > 100:
                continue
                
            comp_str = str(s_comp) + '_' + str(r_comp)
            
            h = open(domain + '_S-R_' + str(s_comp) + '-' + str(r_comp) + '_RESULTS.tsv')
            for i in range(7):
                h.readline()
            header = h.readline()
            for line in h:
                items = line.rstrip().split('\t')
                file_of_origin = items[0]
                if file_of_origin == 'UP000292173_1906665.fasta':   # MANUALLY EXCLUDE MISCLASSIFIED ARCHAEAL PROTEOME (CONFIRMED BY COMMUNICATION WITH UNIPROT)
                    continue
                common_organism_name = items[1]
                junk, uniprot, *junk = items[2].split('|')
                seq = items[3]
                boundaries = items[4]

                if s_comp + r_comp >= threshold:
                    df[file_of_origin] = df.get(file_of_origin, {})  # value is a list, where the first item in the list is another list containing all of the S/R composition combinations in which the protein was identified. The second item in the list is a binary (0 or 1) string value, where 0 indicates that the protein is not an RBP and 1 indicates that the protein is an RBP.
                    df[file_of_origin][uniprot] = df[file_of_origin].get(uniprot, [[], [], [], ''])
                    df[file_of_origin][uniprot][0].append(comp_str)
                    df[file_of_origin][uniprot][1].append(seq)
                    df[file_of_origin][uniprot][2].append(boundaries)
                    df[file_of_origin][uniprot][3] = common_organism_name
                    sr_prots[file_of_origin] = sr_prots.get(file_of_origin, set())
                    sr_prots[file_of_origin].add(uniprot)
            h.close()

    output = open(file_leaders[domain] + '_' + domain + '_SR_proteins_with_Combined_S-R_Above_' + str(threshold) + '.tsv', 'w')
    output.write('\t'.join( ['File of Origin', 'Common Organism Name', 'Protein', 'S-R Composition Thresholds where Protein Identified', 'SR Domain Sequence(s)', 'SR Domain Boundaries'] ) + '\n')
    for file_of_origin in df:
        for uniprot in df[file_of_origin]:
            output.write('\t'.join( [file_of_origin, df[file_of_origin][uniprot][3], uniprot, ','.join(df[file_of_origin][uniprot][0]), ','.join(df[file_of_origin][uniprot][1]), ','.join(df[file_of_origin][uniprot][2])] ) + '\n')
    output.close()
    
    get_full_prot_seqs(sr_prots, domain, threshold)
    

def get_full_prot_seqs(prots, domain, threshold):
    
    seqs_df = {}
    
    for file_of_origin in prots:
        h = open(file_of_origin)
        for seq_record in SeqIO.parse(h, 'fasta'):
            id = str(seq_record.id)
            seq = str(seq_record.seq)
            junk, uniprot, *junk = id.split('|')
            if uniprot in prots[file_of_origin]:
                seqs_df[uniprot] = seq
        h.close()

    pf = open(domain + '_SRproteins_Combined_S-R_Above_' + str(threshold) + '_SequenceDictionary.dat', 'wb')
    pickle.dump(seqs_df, pf)
    pf.close()
    

def get_args(arguments):

    parser = argparse.ArgumentParser(description='Identification of low-complexity domains on the basis of amino acid composition and linear dispersion', prog='LCD-Composer')
    parser.add_argument('domain', help="""Domain of life. Should be "Archaea", "Bacteria", "Eukaryota", or "Viruses".""")
    args = parser.parse_args(arguments)
    
    return args

    
if __name__ == '__main__':
    import sys, argparse, pickle
    from Bio import SeqIO
    taxonid_to_common_name = pickle.load(open('TaxonID_to_CommonOrganismName.dat', 'rb'))
    args = get_args(sys.argv[1:])
    main(args)