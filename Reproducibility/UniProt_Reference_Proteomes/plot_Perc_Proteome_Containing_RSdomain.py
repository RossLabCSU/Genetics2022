
import matplotlib.pyplot as plt
import seaborn as sns

def main():

    plotting_df, human_val = get_data()
    plotting(plotting_df, human_val)

def get_data():

    h = open('PercentageOfProteins_with_RSdomain.tsv')
    header = h.readline()
    
    df = {'Domain':[],
        'Percentage SR prots':[]}
    for line in h:
        domain, file, num_srprots, total_prots, perc_srprots = line.rstrip().split('\t')
        perc_srprots = float(perc_srprots)
        df['Domain'].append(domain)
        df['Percentage SR prots'].append(perc_srprots)
        
        if file == 'UP000005640_9606':
            human_val = perc_srprots
        
    h.close()
    
    return df, human_val
    
    
def plotting(df, human_val):

    colors = sns.color_palette()
    sns.boxplot(x='Domain', y='Percentage SR prots', data=df, palette=['0.8'], showfliers=False, order=['Archaea', 'Bacteria', 'Eukaryota'])
    sns.stripplot(x='Domain', y='Percentage SR prots', data=df, palette=colors, alpha=0.5, order=['Archaea', 'Bacteria', 'Eukaryota'])
        
    plt.xticks(fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel('Domain of Life', fontname='Arial', fontsize=14)
    plt.ylabel('Percentage of Proteins\nwith RS Domain', fontname='Arial', fontsize=14)

    plt.scatter(2, human_val, color='0.2', marker='^', s=35, zorder=3)

    fig = plt.gcf()
    fig.set_size_inches(6, 4)
    plt.savefig('Fig5D_Percentage of Proteins with RS Domain.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    
    
    # PLOT VIRUSES ON A SEPARATE PLOT DUE TO LARGE DIFFERENCE IN Y-AXIS SCALE
    sns.boxplot(x='Domain', y='Percentage SR prots', data=df, palette=['0.8'], showfliers=False, order=['Viruses'])
    sns.stripplot(x='Domain', y='Percentage SR prots', data=df, palette=[colors[3]], alpha=0.5, order=['Viruses'])
        
    plt.xticks(fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel('Domain of Life', fontname='Arial', fontsize=14)
    plt.ylabel('Percentage of Proteins\nwith RS Domain', fontname='Arial', fontsize=14)

    fig = plt.gcf()
    fig.set_size_inches(2.5, 4)
    plt.savefig('Fig5D_Viruses_Percentage of Proteins with RS Domain.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    
    
if __name__ == '__main__':
    main()