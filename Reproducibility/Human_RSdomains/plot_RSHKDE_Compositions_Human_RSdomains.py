
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
            
def main():
    
    df = {'Category':[],
        'Amino Acid':[],
        'Composition':[]}
        
    df = get_data(df)
    single_plot(df)

        
def get_data(df):
    
    h = open('TableS1_Human_RSdomains.tsv')
    header = h.readline()
    for line in h:
        items = line.rstrip().split('\t')
        uniprot = items[0]
        is_rbp = int(items[5])
        is_lcprot = int(items[6])
        merged_domains = ''.join(items[10].split(','))
        comps = items[-20:]
        for i in range(len(amino_acids)):
            if is_rbp == 1 and is_lcprot == 0:
                df['Category'].append( 2 )
            else:
                df['Category'].append( is_rbp )
            df['Amino Acid'].append( amino_acids[i] )
            df['Composition'].append( float(comps[i]) )
    h.close()
    
    return df
        
        
def single_plot(df):

    colors = sns.color_palette('Set1')
    colors = [colors[0], colors[3], colors[1]]
    aaoi = list('RSHKDE')
    df['Amino Acid'] = [aaoi.index(aa) if aa in aaoi else 7 for aa in df['Amino Acid']]
    df = pd.DataFrame.from_dict(df)

    small_df = df[df['Amino Acid'].isin(range(6))]
    
    sns.boxplot(x='Amino Acid', y='Composition', data=small_df, hue='Category', hue_order=[1, 2, 0], palette=['0.8'], showfliers=False)

    ax = plt.gca()
    ax.legend_.remove()
    sns.stripplot(x='Amino Acid', y='Composition', data=small_df, hue='Category', hue_order=[1, 2, 0], palette=colors, size=3, dodge=True, jitter=0.2, alpha=0.7)
    
    plt.xticks([x for x in range(6)], labels=aaoi, fontname='Arial', fontsize=14)
    plt.yticks(fontname='Arial', fontsize=14)
    ax.set_xlabel('')
    ax.set_ylabel('')

    legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[0], label='RBP (Known)', markersize=9),
                        Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[1], label='RBP (New)', markersize=9),
                        Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[2], label='Non-RBP', markersize=9)]
    
    legend = plt.legend(handles=legend_elements, loc='upper right', prop={'family':'Arial', 'size':13}, title='Category', handletextpad=0)
    plt.setp(legend.get_title(), fontsize=14)
    plt.xlabel('Amino Acid', fontname='Arial', fontsize=16)
    plt.ylabel('Percent Composition\nwithin RS Domain', fontname='Arial', fontsize=16)
    plt.tight_layout()
    plt.savefig('Fig 2A - Percent composition of RSHKDE in human RS domains.tiff', bbox_inches ='tight', dpi=600)
    plt.close() 
    
            
if __name__ == '__main__':
    main()