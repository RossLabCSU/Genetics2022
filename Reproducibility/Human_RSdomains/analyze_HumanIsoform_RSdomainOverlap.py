
from Bio import SeqIO
import pickle
origprot_to_isoform = pickle.load(open('Uniprot_OrigProt_to_IsoformList_MappingDict.dat', 'rb'))
isoform_to_origprot = pickle.load(open('Uniprot_Isoform_to_OrigProt_OneToOne_MappingDict.dat', 'rb'))
uniprot_to_common = pickle.load(open('Human_Uniprot_to_GeneName_Conversion.dat', 'rb'))

def main():

    seqs_df = get_prot_seqs()
    domains_df = {}
    orig_SR_prots = get_orig_SR_prots()
    iso_SR_prots = get_iso_SR_prots(orig_SR_prots)
    lc_prots = get_LC_prots()
    merged_domains, merged_boundaries = get_merged_domains(domains_df, seqs_df)
    
    output = open('TableS2_Human_RSdomains_IsoformComparison.tsv', 'w')
    
    h = open('SR_proteins_with_Combined_S-R_Above_70_ALL_ISOFORMS.tsv')
    header = h.readline().rstrip().split('\t')
    output.write('\t'.join([header[0], 'Common Gene Name', 'Isoform of', 'Main isoform is in Long & Caceres Proteins?', 'Original Protein Identified as SR Prot?', 'All Isoforms Identified as SR Prot?', 'Perfect Overlap in Existing RS domains Across SR Isoforms?', 'All Isoforms Identified as SR Prots', 'All Associated Isoforms'] + header[1:] + ['Final Merged Domains', 'Merged Domain Boundaries']) + '\n')
    
    for line in h:
        items = line.rstrip().split('\t')
        prot, comps, seqs, boundaries = items
        common_name = uniprot_to_common[prot]
        if prot in orig_SR_prots:
            orig_prot = prot[:]
            orig_prot_in_SRs = 'Not Applicable'
            isoform_of = 'Not Applicable'
            all_isos = [orig_prot] + origprot_to_isoform[orig_prot]
            is_matched, isoform_in_SRs, all_SR_isos = crosscheck_isoforms(orig_prot, orig_prot, iso_SR_prots, merged_domains, orig_SR_prots)
            
            is_in_lc_prots = '0'
            if orig_prot in lc_prots:
                is_in_lc_prots = '1'
            
            output.write('\t'.join([items[0], common_name, isoform_of, is_in_lc_prots, orig_prot_in_SRs, isoform_in_SRs, is_matched, ','.join(all_SR_isos), ','.join(all_isos)] + items[1:] + [','.join(merged_domains[prot]), ','.join([str(x) for x in merged_boundaries[prot]])]) + '\n')

        else:
            orig_prot = isoform_to_origprot[prot]
            if orig_prot in orig_SR_prots:
                orig_prot_in_SRs = '1'
            else:
                orig_prot_in_SRs = '0'
            all_isos = [orig_prot] + origprot_to_isoform[orig_prot]
            isoform_of = isoform_to_origprot[prot]
            is_matched, isoform_in_SRs, all_SR_isos = crosscheck_isoforms(orig_prot, prot, iso_SR_prots, merged_domains, orig_SR_prots)
            if orig_prot_in_SRs == '0' and isoform_in_SRs == '1':   #CATCH CORNER CASE WHERE ALL ISOFORMS ARE SR PROTEINS BUT THE ORIGINAL PROTEIN IS NOT...CONVERTS THE isoform_in_SRs TO '0'
                isoform_in_SRs = '0'
                
            is_in_lc_prots = '0'
            if orig_prot in lc_prots:
                is_in_lc_prots = '1'
            output.write('\t'.join([items[0], common_name, isoform_of, is_in_lc_prots, orig_prot_in_SRs, isoform_in_SRs, is_matched, ','.join(all_SR_isos), ','.join(all_isos)] + items[1:] + [','.join(merged_domains[prot]), ','.join([str(x) for x in merged_boundaries[prot]])]) + '\n')
            
    h.close()
    output.close()

    
def crosscheck_isoforms(orig_prot, isoform, iso_SR_prots, domains, orig_SR_prots):
    """Function that:
        1) checks whether all isoforms are identified as SR proteins.
            -Value is '1' if True (all isoforms are SR proteins).
            -Value is '0' if False (not all isoforms are SR proteins).
            -Value is 'Not Applicable' if the original protein has no associated isoforms.
        2) checks whether the RS domains of all identified SR proteins perfectly match each other.
            -Value is '1' if all of the identified RS domains perfectly match across all associated isoforms.
            -Value is '0' if the RS domain of at least one isoform differs, or if isoforms have different numbers of RS domains.
        3) collects all isoforms with at least on associated RS domain (later written to output file).
    """
    
    if len(origprot_to_isoform[orig_prot]) == 0:
        is_matched = 'Not Applicable'
        isoform_in_SRs = 'Not Applicable'
        return is_matched, isoform_in_SRs, [orig_prot]
        
    iso_check = []
    match_check = []
    if orig_prot in orig_SR_prots:
        all_SR_isos = [orig_prot]
    else:
        all_SR_isos = []
    for iso in origprot_to_isoform[orig_prot]:
        if iso == isoform:
            isoform_in_SRs = '1'
            is_matched = 'Not Applicable'
            all_SR_isos.append(iso)
            match_check.append(is_matched)
        elif iso in iso_SR_prots:
            isoform_in_SRs = '1'
            is_matched = check_domains(domains, orig_prot, iso)
            all_SR_isos.append(iso)
            match_check.append(is_matched)
        else:
            isoform_in_SRs = '0'
        iso_check.append( isoform_in_SRs )
        
    if '0' in iso_check:
        isoform_in_SRs = '0'
    else:
        isoform_in_SRs = '1'
        
    if '0' in match_check:
        is_matched = '0'
    else:
        is_matched = '1'
        
    return is_matched, isoform_in_SRs, all_SR_isos
    
                
def check_domains(domains, prot, iso):
    """This function is designed to check whether all of the identified RS domains for a set of isoforms are perfectly overlapping.
    This detects insertions/deletions, substitutions, and missing/included RS domains across all isoforms.
    This is only run if the isoform (passed as "iso") is in the iso_SR_prots list.
    """
    
    #FIRST USE ORIGINAL PROTEIN AND CHECK RS DOMAINS AGAINST iso RS DOMAINS.
    match_checks = []
    if prot in domains:
        for domain in domains[prot]:
            domain_checks = []  # domain_checks WILL GET ONE VALUE APPENDED EACH TIME AN RS DOMAIN IS COMPARED. IF domain_checks HAS True ANYWHERE, THAT MEANS THAT THE RS domain FOR prot WAS A PERFECT MATCH FOR AT LEAST PART OF AN iso RS domain (SINCE USING "in" OPERATOR, THIS JUST CHECKS IF THE RS domain IS AT LEAST A PARTIAL MATCH).
            if iso in domains:
                for test_domain in domains[iso]:
                    if domain in test_domain:
                        domain_checks.append(True)
                    else:
                        domain_checks.append(False)
                if True in domain_checks:
                    match_checks.append(True)
                else:
                    match_checks.append(False)
        
    # DO THE REVERSE COMPARISON BY MAKING SURE THAT ALL RS DOMAINS FOR iso ARE CONTAINED IN AT LEAST ONE RS domain FOR prot.
    if iso in domains:
        for domain in domains[iso]:
            domain_checks = []
            if prot in domains:
                for test_domain in domains[prot]:
                    if domain in test_domain:
                        domain_checks.append(True)
                    else:
                        domain_checks.append(False)
                if True in domain_checks:   # AT LEAST ONE True IN DOMAIN_CHECKS MEANS THAT THE RS domain WAS CONTAINED IN AT LEAST ONE RS domain FOR prot.
                    match_checks.append(True)
                else:
                    match_checks.append(False)

    is_matched = '1'
    if False in match_checks:
        is_matched = '0'
        
    return is_matched
    
    
def get_orig_SR_prots():
    
    h = open('SR_proteins_with_Combined_S-R_Above_70.tsv')
    
    header = h.readline()
    prots = []
    for line in h:
        items = line.rstrip().split('\t')
        prots.append( items[0] )
        
    h.close()
    
    return prots
    
def get_iso_SR_prots(orig_SR_prots):
    
    iso_SR_prots = []
    h = open('SR_proteins_with_Combined_S-R_Above_70_ALL_ISOFORMS.tsv')
    header = h.readline()
    for line in h:
        items = line.rstrip().split('\t')
        prot = items[0]
        if prot not in orig_SR_prots:
            iso_SR_prots.append(prot)
            
    h.close()
    
    return iso_SR_prots

    
def get_merged_domains(domains_df, seqs_df):
    
    h = open('SR_proteins_with_Combined_S-R_Above_70_ALL_ISOFORMS.tsv')
    header = h.readline()
    merged_boundaries = {}
    for line in h:
        prot, comps, seqs, boundaries = line.rstrip().split('\t')
        
        comps = comps.split(',')
        seqs = seqs.split(',')
        boundaries = boundaries.split(',')
        
        # GET FORMAL SET OF DOMAIN BOUNDARIES AND SORT IN ASCENDING ORDER FOR LATER MERGING
        bounds = set()
        for boundary in boundaries:
            start, end = boundary[1:-1].split('-')
            bounds.add(  (int(start), int(end))  )
        bounds = sorted( list(bounds) )

        final_bounds = merge_overlapping(bounds)
        merged_boundaries[prot] = final_bounds
        domains_df = get_domains(prot, final_bounds, seqs_df, domains_df)
        
    h.close()
    
    return domains_df, merged_boundaries
    
    
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
    
    
def get_domains(prot, bounds, seqs_df, domains_df):

    for bound in bounds:
        start = bound[0] - 1
        end = bound[1]
        domain = seqs_df[prot][start:end]
        domains_df[prot] = domains_df.get(prot, [])
        domains_df[prot].append(domain)
        
    return domains_df
    
    
def get_prot_seqs():

    h = open('UP000005640_9606.fasta')
    df = {}
    for seq_record in SeqIO.parse(h, 'fasta'):
        junk, id, *junk = str(seq_record.id).split('|')
        seq = str(seq_record.seq)
        df[id] = seq
    h.close()
    
    h = open('UP000005640_9606_additional.fasta')
    for seq_record in SeqIO.parse(h, 'fasta'):
        junk, id, *junk = str(seq_record.id).split('|')
        seq = str(seq_record.seq)
        df[id] = seq
    h.close()
        
    return df
    
    
def get_LC_prots():
    
    prots = []
    h = open('All_SRprots_Long_and_Caceres_2009.txt')
    for line in h:
        prots.append( line.rstrip() )
    h.close()
        
    return prots


if __name__ == '__main__':
    main()