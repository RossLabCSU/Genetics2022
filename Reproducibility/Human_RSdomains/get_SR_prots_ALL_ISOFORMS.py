
def main():

    threshold = 70

    df = {}
    sr_prots = {}
    sr_list = []
    for s_comp in range(20, 105, 5):
        for r_comp in range(20, 105, 5):
            if s_comp + r_comp > 100:
                continue
                
            comp_str = str(s_comp) + '_' + str(r_comp)
            sr_prots[comp_str] = set()
            
            h = open('Human_S-R_' + str(s_comp) + '-' + str(r_comp) + '_RESULTS.tsv')
            sr_prots, sr_list, df = get_prots(h, sr_prots, sr_list, s_comp, r_comp, comp_str, threshold, df)
            
            h.close()
            h = open('Human_S-R_' + str(s_comp) + '-' + str(r_comp) + '_RESULTS_Additional_Isoforms.tsv')
            sr_prots, sr_list, df = get_prots(h, sr_prots, sr_list, s_comp, r_comp, comp_str, threshold, df)
            
            if len(sr_prots[comp_str]) == 0:
                continue
    
    output = open('SR_proteins_with_Combined_S-R_Above_' + str(threshold) + '_ALL_ISOFORMS.tsv', 'w')
    output.write('\t'.join( ['Protein', 'S-R Composition Thresholds where Protein Identified', 'Domain Sequences Identified at Each Corresponding S-R Composition Threshold', 'Domain Boundaries at Each Corresponding S-R Composition Threshold'] ) + '\n')
    for prot in df:
        output.write('\t'.join( [prot, ','.join(df[prot][0]), ','.join(df[prot][1]), ','.join(df[prot][2])] ) + '\n')
    output.close()
    
    
def get_prots(h, sr_prots, sr_list, s_comp, r_comp, comp_str, threshold, df):

    for i in range(7):
        h.readline()
    header = h.readline()
    for line in h:
        items = line.rstrip().split('\t')
        junk, uniprot, *junk = items[0].split('|')
        seq = items[1]
        boundaries = items[2]
        sr_prots[comp_str].add(uniprot)
        sr_list.append(uniprot)
        
        if s_comp + r_comp >= threshold:
            df[uniprot] = df.get(uniprot, [[], [], []])  # value is a list, where the first item in the list is another list containing all of the S/R composition combinations in which the protein was identified. The second item in the list is another list containing each of the domain seqs identified for the corresponding composition criteria. The third item in the list is another list containing the domain boundaries for each of the domains identified at each composition criteria.
            df[uniprot][0].append(comp_str)
            df[uniprot][1].append(seq)
            df[uniprot][2].append(boundaries)
    h.close()
    
    return sr_prots, sr_list, df
    
    
if __name__ == '__main__':
    main()