
from Bio import SeqIO
import pickle
uniprot_to_common = pickle.load(open('Human_Uniprot_to_GeneName_Conversion.dat', 'rb'))
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
mods = ['Phosphorylation', 'Methylation', 'Acetylation', 'Ubiquitination']
mod_dict = pickle.load(open('ActiveDriverDB2021_PTMs_mapped_to_UniprotIDs.dat', 'rb'))
gene_conversion = pickle.load(open('EnsemblGeneID_to_UniProtID_converter.dat', 'rb'))

def main():
    
    seqs_df = get_prot_seqs()
    lc_prots = get_LongCaceres_prots()
    pfam_rrm_accs = get_RRM_pfam_accs()
    gerstberger_rbps = get_gerstberger_rbps()
    goterm_rbps = get_goterm_rbps()
    
    for dataset in ['S-R', 'S-K']:
        if dataset == 'S-R':
            domain_type = 'RS'
        else:
            domain_type = 'KS'
        domains_df = {}
        comps_df = {}
        
        if dataset == 'S-R':
            output = open('TableS1_Human_RSdomains.tsv', 'w')
        else:
            output = open('Human_KSdomains.tsv', 'w')
        
        h = open(dataset[0]+dataset[-1] + '_proteins_with_Combined_' + dataset + '_Above_70_TEST.tsv')
        header = h.readline().rstrip().split('\t')
        new_header = ['UniProt Accession', 'Common Gene Name', 'Contains RRM?', 'Is in Gerstberger et al. RBP dataset?', 'Is annotated with RNA-binding GO term?', 'Is RBP (1 in any of columns 3-5)?', 'Is in Long and Caceres dataset?', dataset + ' Composition Thresholds where Protein Identified', 'Domain Sequences Identified at Each Corresponding ' + dataset + ' Composition Threshold', 'Domain Boundaries at Each Corresponding ' + dataset + ' Composition Threshold', 'Final Merged Domains', 'Merged Domain Boundaries', 'PTMs within ' + domain_type + ' domain', domain_type[1]+domain_type[0]+'/'+domain_type+' Raw Count in SR domain', domain_type[1]+domain_type[0]+'/'+domain_type+' Density in ' + domain_type + ' domain', 'Phosphorylation Count in ' + domain_type + ' domain', 'Methylation Count in ' + domain_type + ' domain', 'Acetylation Count in ' + domain_type + ' domain', 'Ubiquitination Count in ' + domain_type + ' domain', 'Phosphorylation Density in ' + domain_type + ' domain', 'Methylation Density in ' + domain_type + ' domain', 'Acetylation Density in ' + domain_type + ' domain', 'Ubiquitination Density in ' + domain_type + ' domain'] + list(amino_acids)
        output.write( '\t'.join(new_header) + '\n')
        for line in h:
            prot, comp_thresholds, seqs, boundaries = line.rstrip().split('\t')
            
            comp_thresholds = comp_thresholds.split(',')
            seqs = seqs.split(',')
            boundaries = boundaries.split(',')
            
            # GET FORMAL SET OF DOMAIN BOUNDARIES AND SORT IN ASCENDING ORDER FOR LATER MERGING
            bounds = set()
            for boundary in boundaries:
                start, end = boundary[1:-1].split('-')
                bounds.add(  (int(start), int(end))  )
            bounds = sorted( list(bounds) )

            comp_sums = []
            for comp in comp_thresholds:
                s_comp, r_comp = comp.split('_')
                comp_sums.append( int(s_comp) + int(r_comp) )

            final_bounds = merge_overlapping(bounds)
            final_domains = get_domains(prot, final_bounds, seqs_df)
            compositions = calc_compositions(prot, final_domains)
            
            has_rrm, is_gerstberger_rbp, is_goterm_rbp, any_rbp, is_in_lc_prots = run_dataset_crosschecks(prot, lc_prots, gerstberger_rbps, goterm_rbps, pfam_rrm_accs)
            
            mod_str, sr_count, sr_density, modtype_counts, modtype_densities = map_PTMs_and_SRdipeps(prot, final_domains, final_bounds)
            
            output.write('\t'.join( [prot, uniprot_to_common[prot], has_rrm, is_gerstberger_rbp, is_goterm_rbp, any_rbp, is_in_lc_prots, ','.join(comp_thresholds), ','.join(seqs), ','.join(boundaries), ','.join(final_domains), ','.join(['('+str(bound[0])+'-'+str(bound[1])+')' for bound in final_bounds]), mod_str, str(sr_count), str(sr_density)] + [str(x) for x in modtype_counts] + [str(x) for x in modtype_densities] + [str(round(x, 2)) for x in compositions]) + '\n') 

        h.close()
        output.close()

    
def map_PTMs_and_SRdipeps(uniprot, domains, bounds):

    mod_hits = set()
    for bound in bounds:
        start, end = bound
        start_ind = int(start)
        end_ind = int(end)
        mods_hits = get_SRdom_mods_v2(uniprot, start_ind, end_ind, mod_hits)

    sr_count, sr_density = count_SRs(domains)
    
    # SORT MODS BASED ON POSITION
    mod_hits = list(mod_hits)
    if len(mod_hits) > 0:
        residues, positions, mod_types = zip(*mod_hits)
        positions = [int(x) for x in positions]
        positions, residues, mod_types = zip(*sorted(zip(positions, residues, mod_types)))
        positions = [str(x) for x in positions]
        mod_hits = list(zip(residues, positions, mod_types))
        
    mod_str = ';'.join( [str(x).replace("'", '') for x in mod_hits] )

    modtype_counts = []
    modtype_densities = []
    for mod_type in mods:
        count = 0
        for mod in mod_hits:
            if mod[2] == mod_type.lower():
                count += 1
        modtype_counts.append(count)
        modtype_densities.append( count / len(''.join(domains)) )

    # OVERRIDE VALUES IF THERE WERE NO PTMs IN THE IDENTIFIED DOMAINS
    if len(mod_hits) == 0:
        mod_str = 'None'
        modtype_counts = ['nan']*4
        modtype_densities = ['nan']*4
    
    return mod_str, sr_count, sr_density, modtype_counts, modtype_densities
        
        
def count_SRs(domains):
    
    sr_count = 0
    for domain in domains:
        sr_count += domain.count('SR')
        sr_count += domain.count('RS')
        
    sr_density = sr_count / len(''.join(domains))
        
    return sr_count, sr_density
    
    
def get_SRdom_mods_v2(prot, start, end, mod_hits):

    if prot in mod_dict:
        for mod in mod_dict[prot]:
            res = mod[0]
            position = int(mod[1])
            mod_type = mod[2]
            if position >= start and position <= end:
                mod_hits.add( (res, position, mod_type) )
                    
    return mod_hits


def run_dataset_crosschecks(uniprot, lc_prots, gerstberger_rbps, goterm_rbps, pfam_rrm_accs):

    has_rrm = '0'
    if uniprot in pfam_rrm_accs:
        has_rrm = '1'
            
    is_gerstberger_rbp = '0'
    if uniprot in gerstberger_rbps:
        is_gerstberger_rbp = '1'
        
    is_goterm_rbp = '0'
    if uniprot in goterm_rbps:
        is_goterm_rbp = '1'
        
    any_rbp = '0'
    if '1' in [is_gerstberger_rbp, is_goterm_rbp, has_rrm]:
        any_rbp = '1'
        
    is_in_lc_prots = '0'
    if uniprot in lc_prots:
        is_in_lc_prots = '1'
        
        
    return is_gerstberger_rbp, is_goterm_rbp, has_rrm, any_rbp, is_in_lc_prots
    

def calc_compositions(prot, domains):

    combined_domain = ''.join(domains)
    comps = []
    for aa in amino_acids:
        comp = combined_domain.count(aa) / len(combined_domain) * 100
        comps.append(comp)

    return comps
    

def get_domains(prot, bounds, seqs_df):

    domains = []
    for bound in bounds:
        start = bound[0] - 1
        end = bound[1]
        domain = seqs_df[prot][start:end]
        domains.append(domain)

    return domains
        
def merge_overlapping(bounds):

    final_bounds = []
    merged = bounds[0]
    bounds = bounds[1:]
    while len(bounds) > 0:
        bound = bounds[0]
        if bound[0] >= merged[0] and bound[0] <= merged[1]:
            merged = ( merged[0], max(merged[1], bound[1]) )
        else:
            final_bounds.append(merged)
            merged = bounds[0]
        bounds.pop(0)
    final_bounds.append(merged)

    return final_bounds
    
    
def get_prot_seqs():

    h = open('UP000005640_9606.fasta')
    df = {}
    for seq_record in SeqIO.parse(h, 'fasta'):
        junk, id, *junk = str(seq_record.id).split('|')
        seq = str(seq_record.seq)
        df[id] = seq
    h.close()
        
    return df
    
    
def get_LongCaceres_prots():
    
    prots = []
    h = open('All_SRprots_Long_and_Caceres_2009.txt')
    for line in h:
        prots.append( line.rstrip() )
    h.close()
        
    return prots
    
    
def get_RRM_pfam_accs():

    h = open('Pfam_Accessions_for_RRMclan_Only.tsv')
    pfam_rrm_accs = []
    for line in h:
        items = line.rstrip().split('\t')
        pfam_acc = items[0]
        pfam_rrm_accs.append(pfam_acc)
        
    h.close()
    
    return pfam_rrm_accs
    
                
def get_gerstberger_rbps():

    h = open('RBPs_TableS1_from_Gerstberger2014.txt')
    header = h.readline()

    rbps = []
    for line in h:
        items = line.rstrip().split('\t')
        ensembl = items[2]
        
        if ensembl in gene_conversion:
            for prot in gene_conversion[ensembl]:
                rbps.append(prot)
            
    h.close()

    return rbps
    

def get_goterm_rbps():

    h = open('Hsapiens_RNAbinding_Proteins_GOid0003723.txt')
    prots = []
    for line in h:
        prots.append(line.rstrip())
    h.close()
        
    return prots


if __name__ == '__main__':
    main()