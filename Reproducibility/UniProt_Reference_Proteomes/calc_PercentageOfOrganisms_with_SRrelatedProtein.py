
import pickle
import matplotlib.pyplot as plt
import seaborn as sns

domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
table_nums = ['TableS8', 'TableS9', 'TableS10', 'TableS11']
numprots_df = pickle.load(open('TotalNumberOfProteins_per_Organism.dat', 'rb'))

def main():

    output = open('PercentageOfProteins_with_RSdomain.tsv', 'w')
    output.write('\t'.join(['Domain of Life', 'Proteome File', 'Number of Proteins with RS Domain', 'Total Number of Proteins in Proteome', 'Percentage of Proteins with RS Domain']) + '\n')
    
    plotting_df = {'Domain':[],
                'Percentage Organisms with SR prot':[]}
    for i, domain in enumerate(domains):
        table_num = table_nums[i]
        num_srprots_df = get_SRprots(domain, table_num)
        total_organisms = 0
        num_organisms_with_srprot = 0
        for file in numprots_df[domain]:
            num_srprots = num_srprots_df[file]
            total_prots = numprots_df[domain][file]
            perc_prots_with_RSdomain = num_srprots / total_prots * 100

            output.write('\t'.join([domain, file] + [str(x) for x in [num_srprots, total_prots, perc_prots_with_RSdomain]]) + '\n')

            total_organisms += 1
            if num_srprots > 0:
                num_organisms_with_srprot += 1
                
        perc_orgs_with_srprot = num_organisms_with_srprot / total_organisms * 100
        plotting_df['Domain'].append(domain)
        plotting_df['Percentage Organisms with SR prot'].append( perc_orgs_with_srprot )

    output.close()
    
    plotting(plotting_df)


def plotting(df):

    colors = sns.color_palette()
    sns.barplot(x='Domain', y='Percentage Organisms with SR prot', data=df, palette=[colors[0]], order=['Archaea', 'Bacteria', 'Eukaryota', 'Viruses'])
    
    plt.xticks(fontname='Arial', fontsize=12, rotation=90)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel('Domain of Life', fontname='Arial', fontsize=14)
    plt.ylabel('Percentage of Organisms with' + r'$\geq$' + '1\nProtein Containing an RS domain', fontname='Arial', fontsize=14)

    fig = plt.gcf()
    fig.set_size_inches(1.5, 4)
    plt.savefig('Fig5A_Percentage of Organisms with SR protein.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    

def get_SRprots(domain, table_num):
    
    df = {}
    h = open(table_num + '_' + domain + '_SR_proteins_with_Combined_S-R_Above_70_MERGED_DOMAINS_WithPfamAnnotations.tsv')
    header = h.readline()
    
    for line in h:
        file, organism, uniprot, *items = line.rstrip().split('\t')
        file = file.replace('.fasta', '')
        if '_additional' in file:   # SKIP ISOFORMS TO PREVENT POTENTIAL SKEWING OF PERCENTAGES.
            continue
        df[file] = df.get(file, set())
        df[file].add(uniprot)
        
    h.close()
    
    for file in numprots_df[domain]:
        if file in df:
            df[file] = len(df[file])
        else:
            df[file] = 0
            
    return df
    
    
if __name__ == '__main__':
    main()