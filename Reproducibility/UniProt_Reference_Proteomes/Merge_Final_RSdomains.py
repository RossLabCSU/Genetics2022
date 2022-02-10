
import pickle
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']

def main():
    
    for domain in domains:

        seqs_df = pickle.load(open(domain + '_SRproteins_Combined_S-R_Above_70_SequenceDictionary.dat', 'rb'))
        domains_df = {}
        comps_df = {}

        output = open(domain + '_SR_proteins_with_Combined_S-R_Above_70_MERGED_DOMAINS.tsv', 'w')
        
        h = open(domain + '_SR_proteins_with_Combined_S-R_Above_70.tsv')
        header = h.readline().rstrip().split('\t')
        new_header = header + ['Final Merged Domains', 'Merged Domain Boundaries']
        output.write( '\t'.join(new_header) + '\n')
        for line in h:
            file_of_origin, common_organism_name, prot, comps, seqs, boundaries = line.rstrip().split('\t')
            comps = comps.split(',')
            seqs = seqs.split(',')
            boundaries=boundaries.split(',')
            
            # GET FORMAL SET OF DOMAIN BOUNDARIES AND SORT IN ASCENDING ORDER FOR LATER MERGING
            bounds = set()
            for boundary in boundaries:
                start, end = boundary.replace('"', '')[1:-1].split('-')
                bounds.add(  (int(start), int(end))  )
            bounds = sorted( list(bounds) )

            final_bounds = merge_overlapping(bounds)
            domains_df = get_domains(prot, final_bounds, seqs_df, domains_df)
            
            output.write('\t'.join( [file_of_origin, common_organism_name, prot, ','.join(comps), ','.join(seqs), ','.join(boundaries), ','.join(domains_df[prot]), ','.join(['('+str(bound[0])+'-'+str(bound[1])+')' for bound in final_bounds])]) + '\n') 

        
def get_domains(prot, bounds, seqs_df, domains_df):

    for bound in bounds:
        start = bound[0] - 1
        end = bound[1]
        domain = seqs_df[prot][start:end]
        domains_df[prot] = domains_df.get(prot, [])
        domains_df[prot].append(domain)
        
    return domains_df
    
        
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


if __name__ == '__main__':
    main()