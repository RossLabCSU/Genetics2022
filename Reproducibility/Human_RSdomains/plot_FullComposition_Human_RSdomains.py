
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
            
def main():
    
    df = {'Category':[],
        'Amino Acid':[],
        'Composition':[]}
        
    df = get_my_data(df)
    df = pd.DataFrame.from_dict(df)

    plotting(df)
 
 
def get_my_data(df):
    
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
    
        
def plotting(df):

    colors = sns.color_palette('Set1')
    colors = [colors[0], colors[3], colors[1]]
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
        sns.boxplot(x='Amino Acid', y='Composition', data=small_df, hue='Category', hue_order=[1, 2, 0], palette=['0.8'], showfliers=False)
        ax.legend_.remove()
        sns.stripplot(x='Amino Acid', y='Composition', data=small_df, hue='Category', hue_order=[1, 2, 0], palette=colors, size=3, dodge=True, jitter=0.2, alpha=0.7)
        
        plt.xticks(fontname='Arial', fontsize=16)
        plt.yticks(fontname='Arial', fontsize=16)
        ax.set_xlabel('')
        ax.set_ylabel('')

        ax.legend_.remove()
        subplot_index += 1

    fig.text(-0.02, 0.5, 'Percent Composition', va='center', rotation='vertical', fontname='Arial', fontsize=18)
    fig.text(0.5, -0.02, 'Amino Acid', ha='center', fontname='Arial', fontsize=18)
    plt.tight_layout()
    plt.savefig('SR Proteins with combined SR at least 70perc_FullCompAnalysis_RBPs_vs_NonRBPs_LCprotsSplit_COMBINED_RBPs.tiff', bbox_inches ='tight', dpi=600)
    plt.close()   
    
            
if __name__ == '__main__':
    main()