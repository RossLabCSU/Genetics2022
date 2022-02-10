
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def main():

    ptm_types = ['Phosphorylation', 'Methylation', 'Acetylation', 'Ubiquitination']
    h = open('TableS1_Human_RSdomains.tsv')
    header = h.readline()
    ptms_df = {'Category':[],
            'xval':[],
            'yval':[]}
    len_srdomain_df = {'Category':[],
                    'xval':[],
                    'yval':[]}
    combined_density_df = {'Category':[],
                    'xval':[],
                    'yval':[]}
    combined_counts_df = {'Category':[],
                    'xval':[],
                    'yval':[]}
    ptms_density_df = {'Category':[],
                    'xval':[],
                    'yval':[]}
    for line in h:
        items = line.rstrip().split('\t')
        is_rbp = items[5]

        is_lcprot = items[6]
        if is_rbp == '1' and is_lcprot == '1':
            is_rbp = 'RBP (Known)'
        elif is_rbp == '1' and is_lcprot == '0':
            is_rbp = 'RBP (New)'
        else:
            is_rbp = 'Non-RBP'
            
        merged_domains = items[10]
        sr_count = int(items[13])
        sr_per_length = float(items[14])
        combined_density_df['Category'].append(is_rbp)
        combined_density_df['xval'].append(0)
        combined_density_df['yval'].append(sr_per_length)
        combined_counts_df['Category'].append(is_rbp)
        combined_counts_df['xval'].append(0)
        combined_counts_df['yval'].append(sr_count)
        len_srdomain_df['Category'].append(is_rbp)
        len_srdomain_df['xval'].append(0)
        len_srdomain_df['yval'].append(len(merged_domains.replace(',', '')))

        ptm_counts = [int(x) for x in items[15:19]]
        ptm_densities = [float(x) for x in items[19:23]]
        for i in range(4):
            ptms_df['Category'].append(is_rbp)
            ptms_df['xval'].append(i)
            ptms_df['yval'].append(ptm_counts[i])
            ptms_density_df['Category'].append(is_rbp)
            ptms_density_df['xval'].append(i)
            ptms_density_df['yval'].append(ptm_densities[i])
        combined_density_df['Category'].append(is_rbp)
        combined_density_df['xval'].append(1)
        combined_density_df['yval'].append(ptm_densities[0])
        combined_counts_df['Category'].append(is_rbp)
        combined_counts_df['xval'].append(1)
        combined_counts_df['yval'].append(ptm_counts[0])

    h.close()
    
    ptms_df = pd.DataFrame.from_dict(ptms_df)
    ptms_density_df = pd.DataFrame.from_dict(ptms_density_df)
    phos_df = ptms_df[ptms_df['xval'] == 0]
    other_ptms_count_df = ptms_df[ptms_df['xval'] != 0]
    other_ptms_density_df = ptms_density_df[ptms_density_df['xval'] != 0]

    ylabels = ['Total Length of RS Domain(s)', 'Total Number of Sites', 'Density\n(Sites per AA)', 'Total Number of Sites', 'Density\n(Sites per AA)']
    filenames = ['Fig 2B - Length_SRdomains', 'Fig 2C - Combined_SR-RS_Counts', 'Fig 2D - Combined_SR-RS_Phos_Density', 'Fig S5A - Non-phosphorylation PTM Counts', 'Fig S5B - Non-phosphorylation PTM Densities']
    datasets = [len_srdomain_df, combined_counts_df, combined_density_df, other_ptms_count_df, other_ptms_density_df]
    for i in range(len(ylabels)):
        boxplots(datasets[i], ylabels[i], filenames[i])
        
   
def boxplots(df, ylabel, filename):
    
    colors = sns.color_palette('Set1')
    colors = [colors[0], colors[3], colors[1]]
    sns.boxplot(x='xval', y='yval', data=df, hue='Category', showfliers=False, hue_order=['RBP (Known)', 'RBP (New)', 'Non-RBP'], palette=['0.8'])
    sns.stripplot(x='xval', y='yval', data=df, hue='Category', hue_order=['RBP (Known)', 'RBP (New)', 'Non-RBP'], color='black', jitter=True, dodge=True, palette=colors, alpha=0.5)

    ax = plt.gca()
    ax.get_legend().remove()
    
    plt.ylabel(ylabel, fontname='Arial', fontsize=14)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel('')
    
    fig = plt.gcf()
    if filename in ['Fig 2C - Combined_SR-RS_Counts', 'Fig 2D - Combined_SR-RS_Phos_Density']:
        fig.set_size_inches(4.5, 4)
        plt.xticks([0,1], labels=['SR/RS', 'Phosphorylation'], fontname='Arial', fontsize=12)
    elif filename in ['Fig S5A - Non-phosphorylation PTM Counts', 'Fig S5B - Non-phosphorylation PTM Densities']:
        fig.set_size_inches(4.5, 4)
        plt.xticks([0,1,2], labels=['Methylation', 'Acetylation', 'Ubiquitination'], fontname='Arial', fontsize=12)
    else:
        plt.xticks([])
        fig.set_size_inches(2.5, 4)

    plt.savefig(filename + '_Boxplots.tiff', bbox_inches='tight', dpi=600)
    plt.close()

            
if __name__ == '__main__':
    main()