
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
import sys
import argparse
file_leaders = {'Archaea':'TableS8', 'Bacteria':'TableS9', 'Eukaryota':'TableS10', 'Viruses':'TableS11'}
pfam_keywords_df = {'Archaea':['Helicase_C_2', 'DEAD_2'], 'Bacteria':['Helicase_C', 'DEAD'], 'Eukaryota':['Helicase_C', 'DEAD'], 'Viruses':['Corona_nucleoca']}
prot_classes = {'Archaea':'SR-relatedHelicases', 'Bacteria':'SR-relatedHelicases', 'Eukaryota':'SR-relatedHelicases','Viruses':'SR-relatedNucleocapsids'}

def main(args):

    domain_of_life = args.domain_of_life
    domain_locs, prot_to_file = get_domain_locations(pfam_keywords_df[domain_of_life], domain_of_life)
    prot_to_length = get_prot_lengths(prot_to_file)
    
    plotting_df = {'Domain Type':[],
                    'Normalized Distance':[]}
    for prot in domain_locs:
        for domain_type in domain_locs[prot]:
            for domain in domain_locs[prot][domain_type]:
                domain_start = domain[0]
                prot_length = prot_to_length[prot]
                
                normed_start_dist = domain_start / prot_length

                plotting_df['Domain Type'].append(domain_type)
                plotting_df['Normalized Distance'].append(normed_start_dist)
                
    pfam_keywords = pfam_keywords_df[domain_of_life] + ['RS domain']
    prot_class = prot_classes[domain_of_life]
    histplot(plotting_df, pfam_keywords, domain_of_life, prot_class)
    
    
def histplot(df, pfam_keywords, domain, prot_class):

    colors = sns.color_palette()
    for index, kw in enumerate(pfam_keywords):
        vals = [df['Normalized Distance'][i] for i in range(len(df['Normalized Distance'])) if df['Domain Type'][i] == kw]
        sns.distplot(vals, color=colors[index])
        
    plt.xticks(fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel('Normalized Distance from\nProtein N-terminus', fontname='Arial', fontsize=14)
    plt.ylabel('Normalized Frequency', fontname='Arial', fontsize=14)

    legend_elements = [Patch(facecolor=colors[i], edgecolor='w', label=pfam_keywords[i]) for i in range(len(pfam_keywords))]

    ax = plt.gca()
    ax.legend(handles=legend_elements, loc='upper right', handletextpad=0.3, prop={'family':'Arial'})

    fig = plt.gcf()
    fig.set_size_inches(5, 4)
    plt.savefig(domain + '_' + prot_class + '_NormedDistFromStart_Histograms.tiff', bbox_inches='tight', dpi=600)
    plt.close()

    
def get_domain_locations(pfam_keywords, domain_of_life):

    file_leader = file_leaders[domain_of_life]
    h = open(file_leader + '_' + domain_of_life + '_SR_proteins_with_Combined_S-R_Above_70_MERGED_DOMAINS_WithPfamAnnotations.tsv')
    header = h.readline()

    df = {}
    prot_to_file = {}
    for line in h:
        items = line.rstrip().split('\t')
        file = items[0]
        prot = items[2]
        pfams = items[9].split(';')
        if pfam_keywords[0] not in pfams:
            continue
        df[prot] = df.get(prot, {kw:[] for kw in pfam_keywords})
        df[prot]['RS domain'] = df[prot].get('RS domain', [])
        pfam_bounds = items[10].split(';')
        pfam_bounds = [bounds[1:-1].split('-') for bounds in pfam_bounds]
        pfam_bounds = [(int(bounds[0]), int(bounds[1])) for bounds in pfam_bounds]
        
        merged_sr_bounds = items[7].split(',')
        merged_sr_bounds = [bounds[1:-1].split('-') for bounds in merged_sr_bounds]
        merged_sr_bounds = [(int(bounds[0]), int(bounds[1])) for bounds in merged_sr_bounds]
        df[prot]['RS domain'] = merged_sr_bounds
        prot_to_file[prot] = file

        for i, pfam in enumerate(pfams):
            if pfam in pfam_keywords:
                bounds = pfam_bounds[i]
                df[prot][pfam].append(bounds)
                
    h.close()
    
    return df, prot_to_file
    
    
def get_prot_lengths(prot_to_file):
    
    df = {}
    for prot in prot_to_file:
        file = prot_to_file[prot]
        h = open(file)
        for seq_record in SeqIO.parse(h, 'fasta'):
            junk, uniprot, *junk = str(seq_record.id).split('|')
            if uniprot != prot:
                continue
            
            seq = str(seq_record.seq)
            if seq[-1] == '*':
                seq = seq[:-1]
            
            df[prot] = len(seq)
            
        h.close()
        
    return df
    
    
def parse_arguments(arguments) :
    
    parser = argparse.ArgumentParser()
    parser.add_argument('domain_of_life', help = 'input the domain of life (Archaea or Bacteria)')
    args = parser.parse_args(arguments)

    return args

    
if __name__ == '__main__' :
    args = parse_arguments(sys.argv[1:])
    main(args)