
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
total_protcounts = [704232, 28357331, 22511652, 476508]
comp_range = [x for x in range(20, 85, 5)]
import matplotlib.pyplot as plt
import seaborn as sns

def main():

    df = {}
    thresholds = {}
    plotting_df = {'Domain':[],
                'Fraction of SR/SR-related Proteins':[],
                'S+R Composition Bin':[]}
    for domain in domains:
        all_detected_prots = set()
        for s_comp in comp_range:
            for r_comp in comp_range:
                if s_comp + r_comp > 100 or s_comp + r_comp < 70:
                    continue

                df, all_detected_prots = get_hitprots(domain, df, s_comp, r_comp, all_detected_prots)
                
        best_SRbins = [df[prot] for prot in all_detected_prots]
        
        comp_bins = [x for x in range(70, 105, 5)]
        fracs = [best_SRbins.count(x) / len(all_detected_prots) for x in comp_bins]
        
        for i, frac in enumerate(fracs):
            plotting_df['Domain'].append(domain)
            plotting_df['Fraction of SR/SR-related Proteins'].append(frac)
            plotting_df['S+R Composition Bin'].append(comp_bins[i])
                
    barplot(plotting_df)


def barplot(df):
    
    sns.barplot(x='S+R Composition Bin', y='Fraction of SR/SR-related Proteins', data=df, hue='Domain')
    plt.xticks(fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel('Highest S+R Composition Bin', fontname='Arial', fontsize=14)
    plt.ylabel('Fraction of Putative\nSR/SR-related Proteins', fontname='Arial', fontsize=14)
    fig = plt.gcf()
    fig.set_size_inches(5, 4)
    plt.savefig('Fig5B_Fraction_SRprots_PerBin_DomainsOfLife.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    
    
def get_hitprots(domain, df, s_comp, r_comp, all_detected_prots):
    
    h = open(domain + '/' + domain + '_S-R_' + str(s_comp) + '-' + str(r_comp) + '_RESULTS.tsv')
    for i in range(8):
        h.readline()
        
    for line in h:
        items = line.rstrip().split('\t')
        junk, prot, *junk = items[2].split('|')
        all_detected_prots.add(prot)
        df[prot] = df.get(prot, s_comp+r_comp)
        if s_comp + r_comp > df[prot]:
            df[prot] = s_comp+r_comp

    h.close()
    
    return df, all_detected_prots
    
    
if __name__ == '__main__':
    main()