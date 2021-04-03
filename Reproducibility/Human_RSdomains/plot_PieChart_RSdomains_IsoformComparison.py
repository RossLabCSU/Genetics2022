
import matplotlib.pyplot as plt
import seaborn as sns

def main():

    h = open('TableS2_Human_RSdomains_IsoformComparison.tsv')
    header = h.readline()
    match_and_AllIsos = 0
    match_but_NotAllIsos = 0
    notMatch_and_NotAllIsos = 0
    notMatch_but_AllIsos = 0
    one_isoform = 0
    track_prots = set()
    for line in h:
        items = line.rstrip().split('\t')
        common_name = items[1]
        is_origprot_sr = items[4]
        are_allisos_sr = items[5]
        sr_overlap = items[6]
        if common_name not in track_prots:
            if sr_overlap == '1' and are_allisos_sr == '1':
                match_and_AllIsos += 1
            elif sr_overlap == '1' and are_allisos_sr == '0':
                match_but_NotAllIsos += 1
            elif sr_overlap == '0' and are_allisos_sr == '0':
                notMatch_and_NotAllIsos += 1
            elif sr_overlap == '0' and are_allisos_sr == '1':
                notMatch_but_AllIsos += 1
            elif are_allisos_sr == 'Not Applicable':
                one_isoform += 1
        track_prots.add(common_name)

    sizes = [match_and_AllIsos, match_but_NotAllIsos, notMatch_but_AllIsos, notMatch_and_NotAllIsos, one_isoform]
    labels = ['(+) All Isoforms Have RS Domain\n(+) Perfect RS Domain Match', '(-) All Isoforms Have RS Domain\n(+) Perfect RS Domain Match', '(+) All Isoforms Have RS Domain\n(-) Perfect RS Domain Match', '(-) All Isoforms Have RS Domain\n(-) Perfect RS Domain Match', 'Only One Isoform']
    colors = sns.color_palette("Set2")
    total = sum(sizes)
    wedges, texts, autotexts = plt.pie(sizes, labels=labels, autopct=lambda p: '{:.0f}'.format(p * total / 100), startangle=90, textprops={'fontfamily':'Arial', 'fontsize':14}, colors=colors)
    ax = plt.gca()
    ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    legend = ax.legend(wedges, labels,
              title="Category",
              loc="center left",
              bbox_to_anchor=(0.82, 0, 0.5, 1),
              labelspacing=1)

    plt.setp(autotexts, size=14, weight="bold")
    plt.setp(texts, size=1, color='white')
    plt.setp(legend.get_texts(), fontsize=9, va='bottom')
    fig = plt.gcf()
    fig.set_size_inches(6, 4)
    plt.savefig('Fig 3A - IsoformBreakdown_PieChart.tiff', bbox_inches='tight', dpi=600)
    plt.close()

            
if __name__ == '__main__':
    main()