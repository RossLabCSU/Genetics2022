
import matplotlib.pyplot as plt
import seaborn as sns

def main():

    proportion_matching = get_max_proportionMatching()
    proportion_isos = get_num_isos_with_SR()
    df = {'Category':[],
            'Protein':[],
            'Proportion':[]
            }
    match_props = []
    iso_props = []
    prots = []
    for prot in proportion_matching:
        # SKIP PROTEINS FOR WHICH 100% OF ISOFORMS HAVE AN RS domain AND 100% OF RS domainS PERFECTLY MATCH
        if proportion_matching[prot] > 0.99 and proportion_isos[prot] > 0.99:
            continue
        
        match_props.append(proportion_matching[prot])
        iso_props.append(proportion_isos[prot])
        prots.append(prot)
        
    sum_of_values = [3+iso_props[i] if match_props[i] > 0.99 else match_props[i] + iso_props[i] for i in range(len(match_props))]
    sum_of_values, match_props, iso_props, prots = zip(*sorted(zip(sum_of_values, match_props, iso_props, prots), reverse=True))
    
    perfect_match_df = {'Category':[],
            'Protein':[],
            'Proportion':[]
            }
    nonmatch_df = {'Category':[],
            'Protein':[],
            'Proportion':[]
            }
    output = open('ProportionMatching_and_ProportionIsosWithRSdomain_RESULTS.tsv', 'w')
    output.write('\t'.join(['Protein', 'Proportion of RS domains that Match', 'Proportion of Isoforms with RS domain']) + '\n')
    nonmatch_matchprops = []
    nonmatch_isoprops = []
    nonmatch_prots = []
    for i in range(len(match_props)):
        df['Category'].append('RS Domains Align')
        df['Protein'].append(prots[i])
        df['Proportion'].append(match_props[i])
        
        df['Category'].append('Isoforms with RS Domain')
        df['Protein'].append(prots[i])
        df['Proportion'].append(iso_props[i])
        
        if match_props[i] > 0.99:
            perfect_match_df['Category'].append('Isoforms with RS Domain')
            perfect_match_df['Protein'].append(prots[i])
            perfect_match_df['Proportion'].append(iso_props[i])
        else:
            nonmatch_matchprops.append(match_props[i])
            nonmatch_prots.append(prots[i])
            nonmatch_isoprops.append(iso_props[i])

        output.write('\t'.join([str(x) for x in [prots[i], match_props[i], iso_props[i]]]) + '\n')
        
    output.close()

    nonmatch_isoprops, nonmatch_matchprops, nonmatch_prots = zip(*sorted(zip(nonmatch_isoprops, nonmatch_matchprops, nonmatch_prots), reverse=True))
    for i in range(len(nonmatch_matchprops)):
        nonmatch_df['Category'].append('RS Domains Align')
        nonmatch_df['Protein'].append(nonmatch_prots[i])
        nonmatch_df['Proportion'].append(nonmatch_matchprops[i])
        
        nonmatch_df['Category'].append('Isoforms with RS Domain')
        nonmatch_df['Protein'].append(nonmatch_prots[i])
        nonmatch_df['Proportion'].append(nonmatch_isoprops[i])

    plot_nonmatch(nonmatch_df)
    plot_isos_with_perfect_match(perfect_match_df)

    
def plot_nonmatch(df):
    
    colors = sns.color_palette('Set3', 10)
    colors = [colors[0], colors[5]]
    sns.barplot(x='Protein', y='Proportion', hue='Category', hue_order=['Isoforms with RS Domain', 'RS Domains Align'], data=df, palette=colors)
    plt.xticks(fontname='Arial', fontsize=10, rotation=90)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel('Protein', fontname='Arial', fontsize=14)
    plt.ylabel('Proportion', fontname='Arial', fontsize=14)
    plt.ylim(0, 1)
    plt.legend(prop={'size':8})
    fig = plt.gcf()
    fig.set_size_inches(5, 3.5)
    plt.savefig('Fig 3B - ProportionIsoformsWithSRdomain_for_NonMatchProts.tiff', bbox_inches='tight', dpi=600)
    plt.close()

    
def plot_isos_with_perfect_match(df):

    colors = sns.color_palette('Set3', 10)
    colors = [colors[0]]
    sns.barplot(x='Protein', y='Proportion', data=df, palette=colors)
    plt.xticks(fontname='Arial', fontsize=10, rotation=90)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel('Protein', fontname='Arial', fontsize=14)
    plt.ylabel('Proportion of Isoforms\nwith RS Domain', fontname='Arial', fontsize=14)
    plt.ylim(0, 1)
    fig = plt.gcf()
    fig.set_size_inches(13, 3.5)
    plt.savefig('Fig 3C - ProportionIsoformsWithSRdomain_for_IsosWithPerfectSRmatch.tiff', bbox_inches='tight', dpi=600)
    plt.close()

        
def get_max_proportionMatching():
    
    h = open('TableS2_Human_RSdomains_IsoformComparison.tsv')
    header = h.readline()
    seqs_df = {}
    for line in h:
        items = line.rstrip().split('\t')
        prot = items[1]
        domains = items[12]
        seqs_df[prot] = seqs_df.get(prot, [])
        seqs_df[prot].append(domains)
    h.close()
    
    proportion_df = {}
    for prot in seqs_df:
        max_match = 0
        for i in range(len(seqs_df[prot])):
            domain1 = seqs_df[prot][i]
            match = 0.0 - 0.00000001
            for domain2 in seqs_df[prot]:
                if domain1 == domain2:
                    match += 1
            proportion_matching = match / len(seqs_df[prot])
            if proportion_matching > max_match:
                max_match = proportion_matching
        proportion_df[prot] = max_match

    return proportion_df
    
    
def get_num_isos_with_SR():
    
    h = open('TableS2_Human_RSdomains_IsoformComparison.tsv')
    header = h.readline()
    proportion_isos = {}
    for line in h:
        items = line.rstrip().split('\t')
        prot = items[1]
        isos_with_sr = len(items[7].split(','))
        total_isos = len(items[8].split(','))
        proportion = isos_with_sr / total_isos
        domains = items[12]
        proportion_isos[prot] = proportion
    h.close()
    
    return proportion_isos


if __name__ == '__main__':
    main()