
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.lines import Line2D
import math

data_labels = ['Number of significant differentially expressed genes', 'SE', 'MXE', 'A5SS', 'A3SS', 'RI', 'TANDEMUTR']
cats = ['Known', 'New', 'non-SR-related']

def main():

    commons_known, commons_new, ensgs = get_known_srprots()
    
    files = ['SupplementaryData5_DifferentialExpression_and_AlternativeSplicing_K562.txt', 'SupplementaryData5_DifferentialExpression_and_AlternativeSplicing_HepG2.txt']
    labels = ['K562', 'HepG2']
    
    combined_splice_effect_df = {'Splicing Effect':[],
                            'Value':[],
                            'Protein Category':[]}
    for i, file in enumerate(files):
        data = get_data(file)
        label = labels[i]
        df = pd.read_csv(file, sep='\t', skiprows=1)
        df.drop(df.columns[[x for x in range(1,9)]], axis=1, inplace=True)
        df.set_index('RBP', inplace=True)

        new_df = {'Splicing Effect':[],
                'Value':[],
                'Protein Category':[]}
                
        for prot in data:
            if prot in commons_new:
                cat = 'New'
            elif prot in commons_known:
                cat = 'Known'
            else:
                cat = 'non-SR-related'
                
            vals = data[prot]
            for j, val in enumerate(vals):
                if j == 0:
                    continue
                splicing_effect = data_labels[j]
                new_df['Splicing Effect'].append(splicing_effect)
                if val != 0:
                    new_df['Value'].append(math.log10(val))
                else:
                    new_df['Value'].append(-1)
                new_df['Protein Category'].append(cat)
                
            combined_splice_effect_df['Splicing Effect'].append(label)
            combined_splice_effect_df['Value'].append(vals[0])
            combined_splice_effect_df['Protein Category'].append(cat)
                
        plot_splicing_effects(new_df, label, data_labels)

    plot_differentially_expressed(combined_splice_effect_df, 'Cell Line', labels)
    

def plot_differentially_expressed(df, annotation_type, labels):
    
    colors = sns.color_palette('Set1')
    colors = [colors[0], colors[3], colors[8]]
    
    sns.boxplot(x='Splicing Effect', y='Value', data=df, hue='Protein Category', hue_order=cats, palette=['0.8'], showfliers=False, order=labels)
    sns.stripplot(x='Splicing Effect', y='Value', data=df, hue='Protein Category', hue_order=cats, palette=colors, dodge=True, order=labels, alpha=0.5)
    ax = plt.gca()
    
    legend_elements = [Line2D([0], [0], marker='o', color='w', lw=1, label=cats[i], markerfacecolor=colors[i], markersize=7) for i in range(len(cats))]
    ax.legend(handles=legend_elements, loc='upper left', handletextpad=0.0, prop={'size': 9.5}, labelspacing=.2)
    
    if annotation_type == 'Cell Line':
        plt.xticks([i for i in range(len(labels))], labels=labels, fontname='Arial', fontsize=12)
    else:
        plt.xticks([i for i in range(len(labels))], labels=labels, fontname='Arial', fontsize=12, rotation=90)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel(annotation_type, fontname='Arial', fontsize=14)
    plt.ylabel('Number of Differentially\nExpressed Genes', fontname='Arial', fontsize=14)
    plt.ylim(-200, 6900)
    fig = plt.gcf()
    fig.set_size_inches(3.25, 3.5)
    plt.savefig('FigS7D_' + annotation_type + '_DifferentiallyExpressed_Barplot.tif', bbox_inches='tight', dpi=600)
    plt.close()
    

def plot_splicing_effects(df, annotation_type, labels):
    
    colors = sns.color_palette('Set1')
    colors = [colors[0], colors[3], colors[8]]
    
    sns.boxplot(x='Splicing Effect', y='Value', data=df, hue='Protein Category', hue_order=cats, palette=['0.8'], showfliers=False, order=labels[1:])
    sns.stripplot(x='Splicing Effect', y='Value', data=df, hue='Protein Category', hue_order=cats, palette=colors, dodge=True, order=labels[1:], alpha=0.5, size=3)
    ax = plt.gca()

    legend_elements = [Line2D([0], [0], marker='o', color='w', lw=4, label=cats[i], markerfacecolor=colors[i], markersize=8) for i in range(len(cats))]
    ax.legend(handles=legend_elements, loc='upper right', handletextpad=0.0, labelspacing=0.1)

    plt.xticks([i for i in range(len(labels[1:]))], labels=labels[1:], fontname='Arial', fontsize=12, rotation=90)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel(annotation_type, fontname='Arial', fontsize=14)
    plt.ylabel('Log(Number of Genes Affected)', fontname='Arial', fontsize=14)
    plt.ylim(-1.2, 5.9)
    fig = plt.gcf()
    fig.set_size_inches(4, 3.5)
    plt.savefig('FigS7C_' + annotation_type + '_SplicingEffects_Boxplot.tif', bbox_inches='tight', dpi=600)
    plt.close()
    

def get_data(file):

    h = open(file)
    h.readline()
    header = h.readline()
    df = {}
    for line in h:
        items = line.rstrip().split('\t')
        prot = items[0]
        data = [int(x) for x in items[9:]]
        df[prot] = data
    h.close()
    
    return df


def get_known_srprots():

    df = {}
    common_to_ensg = {}
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
            df[uniprot] = 'Known'
            commons_known.add(common)
        else:
            df[uniprot] = 'New'
            commons_new.add(common)
        common_to_ensg[common] = ensg
        
        ensg = ensg.split(', ')
        for name in ensg:
            ensgs.add(name)

    h.close()

    return commons_known, commons_new, ensgs
    
    
if __name__ == '__main__':
    main()
