
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
ensembl_to_uniprot = pickle.load(open('EnsemblGeneID_to_UniProtID_converter.dat', 'rb'))

def main():

    is_known_df, common_to_uniprot = get_known_srprots()
    h = open('SupplementaryInfo3_RNAtarget_Categories.txt')
    header = h.readline()
    known_SRprot_RNAtargets = []
    new_SRprot_RNAtargets = []
    all_other_rbps_RNAtargets = []
    all_targets = set()
    for line in h:
        items = line.rstrip().split('\t')
        ensg = items[2]
        uniprot = [ensg[:]]
        if ensg in ensembl_to_uniprot:
            uniprot = ensembl_to_uniprot[ensg]
        rna_target = items[10]

        all_targets.add(rna_target)
        
        if uniprot[0] in is_known_df:
            if is_known_df[uniprot[0]] == 'Known':
                known_SRprot_RNAtargets.append(rna_target)
            else:
                new_SRprot_RNAtargets.append(rna_target)
        else:
            all_other_rbps_RNAtargets.append(rna_target)
        
    h.close()
    
    df = calculate_frequencies(known_SRprot_RNAtargets, new_SRprot_RNAtargets, all_other_rbps_RNAtargets, all_targets)
    barplot(df)
    
    
def barplot(df):

    colors = sns.color_palette('Set1')
    colors = [colors[0], colors[3], colors[8]]
    
    sns.barplot(x='RNA target', y='Percentage of Proteins', data=df, hue='Category', palette=colors, order=['mRNA', 'snRNA', 'unknown', 'ribosome', 'tRNA', 'rRNA', 'ncRNA', 'diverse', 'snoRNA'])

    plt.xticks(fontname='Arial', fontsize=12, rotation=90)
    plt.yticks([x for x in range(0, 90, 10)], fontname='Arial', fontsize=11)
    plt.xlabel('RNA Target', fontname='Arial', fontsize=14)
    plt.ylabel('Percentage of Proteins', fontname='Arial', fontsize=14)
    fig = plt.gcf()
    fig.set_size_inches(4, 3.5)
    plt.savefig('Fig 4E - RNA target_Barplot_GerstbergerSet_Known_vs_New_vs_AllRBPs.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    
    
def calculate_frequencies(known_SRprot_RNAtargets, new_SRprot_RNAtargets, all_other_rbps_RNAtargets, all_targets):

    df = {'RNA target':[],
        'Percentage of Proteins':[],
        'Category':[]}
        
    ds_labels = ['Known', 'New', 'Non-SR-related']
        
    for target in all_targets:
        for i, dataset in enumerate([known_SRprot_RNAtargets, new_SRprot_RNAtargets, all_other_rbps_RNAtargets]):
            perc = dataset.count(target) / len(dataset) * 100
            df['RNA target'].append(target)
            df['Percentage of Proteins'].append(perc)
            df['Category'].append(ds_labels[i])
            
    return df
    
    
def get_known_srprots():
    
    sr_rbps = get_sr_RBPs()
    
    df = {}
    common_to_uniprot = {}
    file = 'SRrelatedProteins_NameMapping.tsv'
    
    h = open(file)
    header = h.readline()
    for line in h:
        items = line.rstrip().split('\t')
        uniprot = items[0]
        if file.startswith('SRrelatedProteins') and uniprot not in sr_rbps:
            continue
        common = items[1]
        is_known = items[-1]
        if is_known == '1':
            df[uniprot] = 'Known'
        else:
            df[uniprot] = 'New'
        common_to_uniprot[common] = uniprot
    
    h.close()

    return df, common_to_uniprot
    
    
def get_sr_RBPs():

    h = open('TableS1_Human_RSdomains.tsv')
    header = h.readline()
    prots = []
    for line in h:
        items = line.rstrip().split('\t')
        is_rbp = items[5]
        if is_rbp != '1':
            continue
            
        prots.append(items[0])
        
    h.close()
    
    return prots
    
    
if __name__ == '__main__':
    main()