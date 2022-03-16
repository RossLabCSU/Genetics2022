
import pickle
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']

def main():

    residues = ['R', 'S']
    
    for aa in residues: 
        sr_pfam_df, sr_prots_df = get_SRprot_Pfams({domain:[] for domain in domains})
        other_prots_df = get_other_prots(aa)
        
        filtered_sr_prots, filtered_other_prots, shared_prots_df = crosscheck_prots_dfs(sr_prots_df, other_prots_df)
        
        filtered_sr_pfam_df, filtered_sr_prots_df = get_SRprot_Pfams(shared_prots_df)
        filtered_other_pfam_df = get_other_Pfams(aa, shared_prots_df)
        
        sr_top10 = get_top10_pfams(sr_pfam_df)

        # COMPARISON GROUP
        pf = open(aa+'richOnly_proteins_All_Pfams_dict_SHARED_PROTS_EXCLUDED.dat', 'wb')
        pickle.dump(filtered_other_pfam_df, pf)
        pf.close()

        # SR PROTEINS
        pf = open('SR_proteins_top10_Pfams_dict.dat', 'wb')
        pickle.dump(sr_top10, pf)
        pf.close()
        
        pf = open('SR_proteins_All_Pfams_dict_SHARED_PROTS_EXCLUDED.dat', 'wb')
        pickle.dump(filtered_sr_pfam_df, pf)
        pf.close()

    
def get_other_prots(aa):

    prots_df = {}
    for i in range(len(domains)):
        domain = domains[i]
        prots_df[domain] = set()
        h = open(domain + '_' + aa + '-richOnly_WithPfamAnnotations.tsv')
        header = h.readline()
        
        for line in h:
            items = line.rstrip().split('\t')
            junk, uniprot, *junk = items[2].split('|')
            prots_df[domain].add(uniprot)
            
        h.close()
        
    return prots_df
                
                
def crosscheck_prots_dfs(sr_prots_df, other_prots_df):
    
    shared_prots_df = {}
    for domain in domains:
        shared_prots_df[domain] = set()
        for prot in list(sr_prots_df[domain]):
            if prot in other_prots_df[domain]:
                sr_prots_df[domain].remove(prot)
                other_prots_df[domain].remove(prot)
                shared_prots_df[domain].add(prot)
                
    return sr_prots_df, other_prots_df, shared_prots_df
    
                
def get_top10_pfams(pfam_df):

    top10_df = {}
    for domain in domains:
        pfams = list(set(pfam_df[domain]))
        counts = [pfam_df[domain].count(pfam) for pfam in pfams]
        counts, pfams = zip(*sorted(zip(counts, pfams), reverse=True))
        top10_df[domain] = pfams[:10]
        
    return top10_df
            

def get_other_Pfams(aa, shared_prots_df):
    
    pfam_df = {}
    prots_df = {}
    for i in range(len(domains)):
        domain = domains[i]
        pfam_df[domain] = []
        prots_df[domain] = set()
        h = open(domain + '_' + aa + '-richOnly_WithPfamAnnotations.tsv')
        header = h.readline()
        
        for line in h:
            items = line.rstrip().split('\t')
            junk, uniprot, *junk = items[2].split('|')
            if uniprot in shared_prots_df[domain]:
                continue
            if uniprot not in prots_df[domain]:     # PREVENT DOUBLE-COUNTING FOR PROTEINS WITH MULTIPLE LCDs (EACH LCD IS REPRESENTED ON A SEPARATE ROW)
                prots_df[domain].add(uniprot)
                if items[8] != 'N/A':
                    pfams = items[8].split(';')
                    pfam_df[domain] += list(set(pfams))     # PREVENT EXTRA COUNTING FOR PROTEINS WITH MULTIPLE COPIES OF THE SAME PFAM DOMAIN (e.g. PROTEINS WITH MULTIPLE RRMs). THIS WOULD OTHERWISE SKEW STATISTICS FOR LONG PROTEINS WITH REPEATED DOMAINS.
            
        h.close()
        
    return pfam_df
            
    
def get_SRprot_Pfams(shared_prots_df):
    
    leaders = {'Archaea':'TableS8', 'Bacteria':'TableS9', 'Eukaryota':'TableS10', 'Viruses':'TableS11'}
    pfam_df = {}
    prots_df = {}
    
    for i in range(len(domains)):
        domain = domains[i]
        leader = leaders[domain]
        pfam_df[domain] = []
        prots_df[domain] = set()
        h = open(leader + '_' + domain + '_SR_proteins_with_Combined_S-R_Above_70_MERGED_DOMAINS_WithPfamAnnotations.tsv')
        header = h.readline()
        for line in h:
            items = line.rstrip().split('\t')
            uniprot = items[2]
            if uniprot in shared_prots_df[domain]:
                continue
            prots_df[domain].add(uniprot)
            if items[9] != 'N/A':
                pfams = items[9].split(';')
                pfam_df[domain] += list(set(pfams))
            
        h.close()
        
    return pfam_df, prots_df


if __name__ == '__main__':
    main()
