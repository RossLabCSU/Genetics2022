
import datetime
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D


def main():

    known, new, nonSR = get_combined_prots()
    df = {}
    labels = ['New', 'Known', 'Non-SR-related']
    prots_df = {cat:set() for cat in labels}
    plotting_df = {'Number of Transcripts Bound':[],
                'Category':[]}

    file = 'catrapid_human_basic_normalised_AllProteins_SCORE_ATLEAST_21p297.txt'
    num_lines = get_filelength(file)
    h = open(file)
    h.readline()

    for line in tqdm(h, total=num_lines):
        prot, transcript, score = line.rstrip().split('\t')
        if prot in new:
            cat = 'New'
        elif prot in known:
            cat = 'Known'
        else:
            cat = 'Non-SR-related'
        df[prot] = df.get(prot, set())
        df[prot].add(transcript)
        prots_df[cat].add(prot)

    h.close()

    for label in labels:
        plotting_df = calc_num_transcripts_bound(df, prots_df, label, plotting_df)

    print('Largest number of transcripts bound:', max(plotting_df['Number of Transcripts Bound']))
    print('Performing plotting...')
    
    plot_numtranscriptsbound(plotting_df)
    
                
def plot_numtranscriptsbound(plotting_df):
    colors = sns.color_palette('Set1')
    colors = [colors[0], colors[3], colors[8]]
    
    temp_df = plotting_df.copy()
    temp_df['Number of Transcripts Bound'] = [plotting_df['Number of Transcripts Bound'][i] if plotting_df['Category'][i] != 'Non-SR-related' else -20000 for i in range(len(plotting_df['Number of Transcripts Bound']))]
    sns.boxplot(x='Category', y='Number of Transcripts Bound', data=temp_df, showfliers=False, palette=['0.8'], order=['Known', 'New', 'Non-SR-related'])
    sns.stripplot(x='Category', y='Number of Transcripts Bound', data=temp_df, palette=colors, alpha=0.5, size=4, order=['Known', 'New', 'Non-SR-related'])
    
    temp_df = plotting_df.copy()
    temp_df['Number of Transcripts Bound'] = [plotting_df['Number of Transcripts Bound'][i] if plotting_df['Category'][i] == 'Non-SR-related' else -20000 for i in range(len(plotting_df['Number of Transcripts Bound']))]
    sns.boxplot(x='Category', y='Number of Transcripts Bound', data=temp_df, showfliers=False, palette=['0.8'], order=['Known', 'New', 'Non-SR-related'])
    sns.stripplot(x='Category', y='Number of Transcripts Bound', data=temp_df, palette=colors, alpha=0.2, size=2, order=['Known', 'New', 'Non-SR-related'])

    labels = ['Known', 'New', 'non-SR-related']
    ax = plt.gca()
    legend_elements = [Line2D([0], [0], marker='o', color='w', lw=1, label=labels[i], markerfacecolor=colors[i], markersize=7) for i in range(len(labels))]
    ax.legend(handles=legend_elements, loc='upper left', handletextpad=0.0, prop={'size': 9.5}, labelspacing=.2)
    
    plt.xticks([i for i in range(len(labels))], labels=labels, fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel('Protein Category', fontname='Arial', fontsize=14)
    plt.ylabel('Predicted Number of\nTranscripts Bound', fontname='Arial', fontsize=14)
    plt.ylim(-2000, 59000)
    
    fig = plt.gcf()
    fig.set_size_inches(4.5, 3.5)
        
    plt.savefig('Fig 4H - NumberOfTranscriptsBound.tif', bbox_inches='tight', dpi=600)
    plt.close()
                
                
def calc_num_transcripts_bound(df, prots_df, label, plotting_df):

    for prot in prots_df[label]:
        transcripts = df[prot]
        t_bound_A = len(transcripts)
        plotting_df['Number of Transcripts Bound'].append(t_bound_A)
        plotting_df['Category'].append(label)

    return plotting_df
                

def get_filelength(file):
    
    h = open(file)
    count = 0
    for line in h:
        count += 1
        
    h.close()
    
    return count
    
def get_combined_prots():

    h = open('Combined_SRprots_and_mRNA-bindingProts_UniProt_Accessions.tsv')
    header = h.readline()

    known = []
    new = []
    nonSR = []
    for line in h:
        prot, cat = line.rstrip().split('\t')
        if cat == 'non-SR-related':
            nonSR.append(prot)
        elif cat == 'New SR-related':
            new.append(prot)
        else:
            known.append(prot)
            
    h.close()

    return known, new, nonSR


if __name__ == '__main__':
    main()