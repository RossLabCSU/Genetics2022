
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import statistics

function_labels = ['Essential Genes', 'Splicing regulation', 'Spliceosome', 'RNA modification', '3\' end processing', 'rRNA processing', 'Ribosome & basic transl.', 'RNA stability & decay', 'microRNA processing', 'RNA localization', 'RNA export', 'Translation regulation', 'tRNA regulation', 'mtRNA regulation', 'Viral RNA regulation', 'snoRNA/snRNA/telomerase', 'P-body/stress granules', 'Exon Junction Complex', 'Novel RBP', 'Other']
loc_labels = ['Nuclei', 'Nucleolus', 'Speckles', 'PML bodies', 'Cajal bodies', 'Cytoplasm', 'Mitochondria', 'Golgi', 'P bodies', 'ER', 'Cytoskeleton', 'Microtubule', 'Actin', 'Nuclear release mitosis', 'Cell Cortex']
pfam_labels = ['RRM', 'ZNF', 'KH', 'Helicase', 'Nuclease', 'dRBM', 'PUM_HD']
experiment_labels = ['eCLIP HepG2', 'eCLIP K562', 'RNAseq HepG2', 'RNAseq K562', 'RBNS', 'Image HepG2', 'ChIPseq HepG2', 'ChIPseq K562']
cats = ['Known', 'New', 'non-SR-related']

all_labels = [function_labels, loc_labels, pfam_labels, experiment_labels]
annotation_types = ['Functions', 'Localizations', 'Pfams', 'Experiments']
file_leaders = ['Fig4F_', 'Fig4C_', 'Fig4D_', 'Fig4A_']

def main():

    commons_known, commons_new, ensgs = get_known_srprots()
    commons = commons_new & commons_known

    df = pd.read_csv('SupplementaryData1_ENCORE_data.txt', sep='\t')
    df.drop(df.columns[[x for x in range(22, 52)]], axis=1, inplace=True)
    df.drop('GeneID', axis=1, inplace=True)
    df = df.replace(np.nan, -1)
    df.set_index('RBP name', inplace=True)

    data = get_data()

    plotting_df = {'Annotation Type':[],    
                'Protein Category':[],
                'Average':[]}
    for key in data:
        if key == 'Functions':  # THIS SPECIFICALLY EXCLUDES A FEW OF THE COLUMNS THAT DO NOT REPRESENT SPECIFIC FUNCTIONS (namely, "Essential Genes", "Novel RBP", and "Other")
            known = [data[key][prot][1:-2] for prot in data[key] if prot in commons_known]
            new = [data[key][prot][1:-2] for prot in data[key] if prot in commons_new]
            non_sr = [data[key][prot][1:-2] for prot in data[key] if prot not in commons_known and prot not in commons_new]
        else:
            known = [data[key][prot] for prot in data[key] if prot in commons_known]
            new = [data[key][prot] for prot in data[key] if prot in commons_new]
            non_sr = [data[key][prot] for prot in data[key] if prot not in commons_known and prot not in commons_new]
        known = [[int(x) for x in l if x != '#N/A'] for l in known]
        known_ave = statistics.mean([sum(x) for x in known])
        
        plotting_df['Annotation Type'].append(key)
        plotting_df['Protein Category'].append('Known')
        plotting_df['Average'].append(known_ave)

        new = [[int(x) for x in l if x != '#N/A'] for l in new]
        new_ave = statistics.mean([sum(x) for x in new])
        
        plotting_df['Annotation Type'].append(key)
        plotting_df['Protein Category'].append('New')
        plotting_df['Average'].append(new_ave)

        non_sr = [[int(x) for x in l if x != '#N/A'] for l in non_sr]
        non_sr_ave = statistics.mean([sum(x) for x in non_sr])
        
        plotting_df['Annotation Type'].append(key)
        plotting_df['Protein Category'].append('non-SR-related')
        plotting_df['Average'].append(non_sr_ave)

    plot_mean_num_annotations(plotting_df)
    
    new_df = {'Protein Category':[],
            'Value':[],
            'Localization Category':[]}
            
    new_df = {'Known':[], 'New':[], 'non-SR-related':[]}
    prots_df = {'Known':[], 'New':[], 'non-SR-related':[]}
    for prot in data['Localizations']:
        locs = [loc_labels[i] for i in range(len(data['Localizations'][prot])) if data['Localizations'][prot][i] == '1']
        if prot in commons_new:
            new_df['New'] += locs
            prots_df['New'].append(prot)
        elif prot in commons_known:
            new_df['Known'] += locs
            prots_df['Known'].append(prot)
        else:
            new_df['non-SR-related'] += locs
            prots_df['non-SR-related'].append(prot)

    for i, labels in enumerate(all_labels):
        annotation_type = annotation_types[i]
        file_leader = file_leaders[i]

        new_df = {'Category':[],
                'Percentage of Proteins':[],
                'Protein Category':[]}
                
        for j, label in enumerate(labels):
            temp = {cat:0 for cat in cats}
            prots_df = {'Known':[], 'New':[], 'non-SR-related':[]}
            for prot in data[annotation_type]:
                if prot in commons_new:
                    cat = 'New'
                    prots_df['New'].append(prot)
                elif prot in commons_known:
                    cat = 'Known'
                    prots_df['Known'].append(prot)
                else:
                    cat = 'non-SR-related'
                    prots_df['non-SR-related'].append(prot)
                    
                indicator = data[annotation_type][prot][j]
                if indicator == '1':
                    temp[cat] += 1
            
            percs = {cat:temp[cat] / len(prots_df[cat]) * 100 for cat in temp}

            for cat in percs:
                new_df['Category'].append(label)
                new_df['Percentage of Proteins'].append(percs[cat])
                new_df['Protein Category'].append(cat)

        barplot(new_df, annotation_type, file_leader, labels)


def barplot(df, annotation_type, file_leader, labels):

    colors = sns.color_palette('Set1')
    colors = [colors[0], colors[3], colors[8]]

    known_vals = [df['Percentage of Proteins'][i] for i, cat in enumerate(df['Protein Category']) if cat=='Known']
    known_vals, labels = zip(*sorted(zip(known_vals, labels), reverse=True))
    sns.barplot(x='Category', y='Percentage of Proteins', data=df, hue='Protein Category', palette=colors, order=labels)
    ax = plt.gca()
    ax.legend(loc='upper right')
    plt.xticks([i for i in range(len(labels))], labels=[x.replace('Nuclear release mitosis', 'Nuclear release,\nmitosis') for x in labels], fontname='Arial', fontsize=12, rotation=90, linespacing=0.8)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel(annotation_type[:-1] + ' Category', fontname='Arial', fontsize=14)
    plt.ylabel('Percentage of Proteins', fontname='Arial', fontsize=14)
    fig = plt.gcf()
    if annotation_type in ['Functions', 'Localizations']:
        fig.set_size_inches(5, 3.5)
    else:
        fig.set_size_inches(4.5, 3.5)
    plt.savefig(file_leader + annotation_type + '_Barplot.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    

def plot_mean_num_annotations(df):

    colors = sns.color_palette('Set1')
    colors = [colors[0], colors[3], colors[8]]
    
    sns.barplot(x='Annotation Type', y='Average', data=df, hue='Protein Category', palette=colors)
    plt.xticks(fontname='Arial', fontsize=12, rotation=90)
    plt.yticks(fontname='Arial', fontsize=11)
    plt.xlabel('Annotation Category', fontname='Arial', fontsize=14)
    plt.ylabel('Mean # of Literature-Based\nAnnotations Per Protein', fontname='Arial', fontsize=12)
    ax = plt.gca()
    ax.legend(loc='upper left', handletextpad=0.3, prop={'size':9.5}, labelspacing=0.2)
    plt.ylim(0, 4.5)
    fig = plt.gcf()
    fig.set_size_inches(4, 4)
    plt.savefig('Fig4B_MeanAnnotationNumber.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    
    
def get_data():

    h = open('SupplementaryData1_ENCORE_data.txt')
    h.readline()
    header = h.readline()
    count = 0
    df = {'Functions':{},
            'Localizations':{},
            'Pfams':{},
            'Experiments':{}}
    for line in h:
        items = line.rstrip().split('\t')
        common = items[0]
        functions = items[2:22]
        localization = items[22:37]
        pfam = items[37:44]
        experiments = items[44:]
        
        df['Functions'][common] = df['Functions'].get(common, functions)
        if len(set(localization)) == 1 and '#N/A' in localization:
            pass
        else:
            df['Localizations'][common] = df['Localizations'].get(common, localization)
        df['Pfams'][common] = df['Pfams'].get(common, pfam)
        df['Experiments'][common] = df['Experiments'].get(common, experiments)
        
    h.close()
    
    return df
 
    
def get_known_srprots():

    commons_known = set()
    commons_new = set()
    ensgs = set()
    file = 'SRrelatedProteins_NameMapping.tsv'
    h = open(file)
    header = h.readline()
    for line in h:
        items = line.rstrip().split('\t')
        uniprot = items[0]
        common = items[1]
        ensg = items[2]
        is_known = items[-1]
        if is_known == '1':
            commons_known.add(common)
        else:
            commons_new.add(common)
        
        ensg = ensg.split(', ')
        for name in ensg:
            ensgs.add(name)

    h.close()

    return commons_known, commons_new, ensgs


if __name__ == '__main__':
    main()