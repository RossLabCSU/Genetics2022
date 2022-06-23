
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
file_leaders = {'Archaea':'TableS8', 'Bacteria':'TableS9', 'Eukaryota':'TableS10', 'Viruses':'TableS11'}
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    df = {'Amino Acid':[],
        'Percent Composition':[],
        'Domain of Life':[]}
        
    grouped_df = {'Amino Acid Group':[],
        'Percent Composition':[],
        'Domain of Life':[]}
    
    for domain in domains:
        filename = file_leaders[domain] + '_' + domain + '_SR_proteins_with_Combined_S-R_Above_70_MERGED_DOMAINS_WithPfamAnnotations'
        h = open(filename + '.tsv')
        header = h.readline().rstrip().split('\t')

        for line in h:
            items = line.rstrip().split('\t')
            merged_domains = ''.join(items[6].split(','))
            for aa in amino_acids:
                perc_comp = merged_domains.count(aa) / len(merged_domains) * 100
                df['Amino Acid'].append(aa)
                df['Percent Composition'].append(perc_comp)
                df['Domain of Life'].append(domain)
        h.close()
        
    plotting(df)
    
    
def plotting(df):

    df = pd.DataFrame.from_dict(df)
    colors = sns.color_palette()
    subplot_index = 1
    fig, ax = plt.subplots(4,1,sharex=True, sharey=False, figsize=(8,10))
    big_subplot = fig.add_subplot(111)
    big_subplot.spines['top'].set_color('none')
    big_subplot.spines['bottom'].set_color('none')
    big_subplot.spines['left'].set_color('none')
    big_subplot.spines['right'].set_color('none')
    big_subplot.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    
    for aa_range in [(0,5), (5,10), (10,15), (15,20)]:
        small_df = df[df['Amino Acid'].isin( list(amino_acids[aa_range[0]:aa_range[1]]) )]

        ax = plt.subplot(4,1, subplot_index)

        sns.boxplot(x='Amino Acid', y='Percent Composition', data=small_df, hue='Domain of Life', palette=colors, showfliers=False)

        ax.get_legend().remove()
        
        plt.xticks(fontname='Arial', fontsize=16)
        plt.yticks([x for x in range(0, 80, 10)], fontname='Arial', fontsize=16)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_ylim(-3, 78)

        subplot_index += 1

    fig.text(-0.02, 0.5, 'Percent Composition', va='center', rotation='vertical', fontname='Arial', fontsize=18)
    fig.text(0.5, -0.02, 'Amino Acid', ha='center', fontname='Arial', fontsize=18)
    plt.tight_layout()
    plt.savefig('Fig S7 - Full Composition Analysis of RS Domains - Domains of Life_Boxplot.tiff', bbox_inches ='tight', dpi=600)
    plt.close() 
        

if __name__ == '__main__':
    main()