
import matplotlib.pyplot as plt
import seaborn as sns
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
file_leaders = {'Archaea':'TableS8', 'Bacteria':'TableS9', 'Eukaryota':'TableS10', 'Viruses':'TableS11'}

def main():

    for domain in domains:  
        pfams = get_pfam_annots(domain)
        freqs, pfams = calc_freqs(pfams)
        
        piechart(freqs, pfams, domain)
        

def piechart(freqs, pfams, domain):

    # DEFAULT COLOR PALETTE USED FOR EUKARYOTES
    colors = sns.color_palette("Set2")
    colors += sns.color_palette('pastel')
    
    # MANUAL ALTERATION OF COLOR PALETTES TO ENSURE THAT PFAMS APPEARING IN PIE CHARTS FOR MULTIPLE ORGANISMS ARE ASSIGNED A CONSISTENT COLOR
    if domain == 'Archaea':
        color1, color2 = colors[2], colors[0]
        colors[2], colors[0] = color2, color1
        
        color1, color2 = colors[1], colors[4]
        colors[1], colors[4] = color2, color1
    
    elif domain == 'Bacteria':
        color1, color2 = colors[0], colors[1]
        colors[0], colors[1] = color2, color1
        
        color1, color2 = colors[1], colors[2]
        colors[1], colors[2] = color2, color1
        
        color1, color2 = colors[4], colors[2]
        colors[4], colors[2] = color2, color1

    total = sum(freqs)
    other = sum(freqs[10:])
    labels = list(pfams[:10]) + ['Other']
    sizes = list(freqs[:10]) + [other]
    
    if domain in ['Bacteria', 'Eukaryota']:
        wedges, texts, autotexts = plt.pie(sizes, labels=labels, pctdistance=1.12, autopct=lambda p: '{:.0f}'.format(p * total / 100), startangle=90, textprops={'fontfamily':'Arial', 'fontsize':12}, colors=colors)
    else:
        wedges, texts, autotexts = plt.pie(sizes, labels=labels, pctdistance=1.1, autopct=lambda p: '{:.0f}'.format(p * total / 100), startangle=90, textprops={'fontfamily':'Arial', 'fontsize':12}, colors=colors)
    ax = plt.gca()
    ax.axis('equal')  # EQUAL ASPECT RATIO ENSURES THAT PIE IS DRAWN AS A CIRCLE

    legend = ax.legend(wedges, labels,
              title="Pfam Name",
              loc="center left",
              bbox_to_anchor=(0.86, 0, 0.5, 1),
              labelspacing=1)

    leaders = {'Archaea':'Fig4C', 'Bacteria':'Fig4D', 'Eukaryota':'Fig4B', 'Viruses':'Fig4E'}
    plt.setp(texts, size=1, color='white')
    plt.setp(legend.get_texts(), fontsize=8, va='bottom')
    fig = plt.gcf()
    fig.set_size_inches(6, 4)
    plt.savefig(leaders[domain] + '_' + domain + '_Top10FrequentPfams_PieChart.tiff', bbox_inches='tight', dpi=600)
    plt.close()
        
        
def get_pfam_annots(domain):

    pfams = []
    h = open(file_leaders[domain] + '_' + domain + '_SR_proteins_with_Combined_S-R_Above_70_MERGED_DOMAINS_WithPfamAnnotations.tsv')
    header = h.readline()
    
    for line in h:
        items = line.rstrip().split('\t')
        if items[9] == 'N/A':   # SKIP PROTEINS FOR WHICH NO PFAM DOMAINS WERE IDENTIFIED
            continue
        annots = items[9].split(';')
        annots = list(set(annots))  # MAKE ANNOTS NON-REDUNDANT...COUNTS EACH PFAM ONLY ONCE FOR EACH PROTEIN REGARDLESS OF HOW MANY TIMES IT OCCURS IN THAT PROTEIN
        pfams += annots
    
    h.close()
    
    return pfams
    
    
def calc_freqs(pfams):

    annots = []
    freqs = []
    for pfam in set(pfams):
        freq = pfams.count(pfam)
        annots.append(pfam)
        freqs.append(freq)
        
    freqs, annots = zip(*sorted(zip(freqs, annots), reverse=True))

    return freqs, annots
    

if __name__ == '__main__':
    main()
