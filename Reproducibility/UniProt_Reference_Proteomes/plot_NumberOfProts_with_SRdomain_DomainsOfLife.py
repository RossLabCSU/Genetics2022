
import matplotlib.pyplot as plt
import pickle
import seaborn as sns
import pandas as pd
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
taxonid_to_common = pickle.load(open('TaxonID_to_CommonOrganismName.dat', 'rb'))

def main():

    file_df = pickle.load(open('DomainsOfLife_AllOrganism_Filenames.dat', 'rb'))
    file_df['Humans'] = ['UP000005640_9606']
    sr_count_df = get_numprots_with_SRs()
    df = {'Number of SR Proteins':[],
        'Proteome':[]
        }
    output = open('Number_of_SRproteins_Per_Organism.tsv', 'w')
    output.write('\t'.join(['Domain of Life', 'File of Origin','Common Organism Name', 'Number of SR/SR-related Proteins (isoforms included)']) + '\n')
        
    index = 0
    for domain in domains:
        nums_srs = []
        for file in file_df[domain]:
            if 'UP000005640_9606' in file:  # SKIP HUMAN PROTEOME BECAUSE THIS IS LATER PLOTTED MANUALLY WITH A DIFFERENT MARKER
                continue
            if file in sr_count_df[domain]:
                count = sr_count_df[domain][file]
            else:
                count = 0
                
            nums_srs.append(count)
            df['Number of SR Proteins'].append(count)
            df['Proteome'].append(index)

            junk, taxonid = file.replace('.fasta', '').split('_')
            if taxonid in taxonid_to_common:
                common_name = taxonid_to_common[taxonid]
            else:
                common_name = file
            output.write('\t'.join( [domain, file, common_name, str(count)] ) + '\n')

        index += 1
        
    output.close()
    
    human_count = sr_count_df['Eukaryota']['UP000005640_9606']

    plotting(df)
    plot_eukaryota(df, human_count)
    
    
def plot_eukaryota(df, human_count):

    df = pd.DataFrame.from_dict(df)
    small_df = df[df['Proteome'] == 2]
    colors = sns.color_palette()
    
    sns.boxplot(x='Proteome', y='Number of SR Proteins', data=small_df, palette=['0.8'], showfliers=False)
    sns.stripplot(x='Proteome', y='Number of SR Proteins', data=small_df, palette=[colors[2]], alpha=0.5)
    plt.scatter(0, human_count, color='0.2', marker='^', s=35, zorder=3)
    plt.xticks([0], labels=['Eukaryota'], fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel('Domain of Life', fontname='Arial', fontsize=14)
    plt.ylabel('Number of Unique Genes Encoding\nSR/SR-related Proteins', fontname='Arial', fontsize=14)
    fig = plt.gcf()
    fig.set_size_inches(2.5, 4)
    plt.savefig('Fig 4A - Eukaryotes_NumberSRproteins.tiff', bbox_inches='tight', dpi=600)
    plt.close()

    
def plotting(df):

    df = pd.DataFrame.from_dict(df)
    
    small_df = df[df['Proteome'].isin([0,1,3])]
    
    colors = sns.color_palette()
    colors = [colors[i] for i in [0,1,3]]
    sns.boxplot(x='Proteome', y='Number of SR Proteins', data = small_df, palette=['0.8'], showfliers=False)
    sns.stripplot(x='Proteome', y='Number of SR Proteins', data=small_df, palette=colors, alpha=0.5)
    
    plt.xticks([0,1,2], labels=['Archaea', 'Bacteria', 'Viruses'], fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel('Domain of Life', fontname='Arial', fontsize=14)
    # plt.xlabel('')
    plt.ylabel('Number of Unique Genes Encoding\nSR/SR-related Proteins', fontname='Arial', fontsize=14)
    fig = plt.gcf()
    fig.set_size_inches(6, 4)
    
    plt.savefig('Fig 4A - Non-Eukaryotes_NumberSRproteins.tiff', bbox_inches='tight', dpi=600)
    plt.close()

    
def get_numprots_with_SRs():
    
    df = {}
    for domain in domains:
        df[domain] = {}
        h = open(domain + '_SR_proteins_with_Combined_S-R_Above_70.tsv')
        header = h.readline()
        for line in h:
            items = line.rstrip().split('\t')
            file_of_origin, ext = items[0].split('.')
            if 'additional' in file_of_origin:
                temp = file_of_origin.split('_')
                file_of_origin = temp[0] + temp[1]
            df[domain][file_of_origin] = df[domain].get(file_of_origin, 0) + 1
            
        h.close()

    return df
    
if __name__ == '__main__':
    main()