
import matplotlib.pyplot as plt

def main():
    
    # HARD-CODED VALUES TAKEN FROM FIG S2A.
    # NOTE: FOR BOTH 80% AND 85% COMBINED S+R COMPOSITION THRESHOLDS, THERE IS ONE PROTEIN THAT IS NOT DEPICTED DUE TO GEOMETRIC CONSTRAINTS WHEN DRAWING THE VENN DIAGRAMS: THIS PROTEIN IS IDENTIFIED IN THIS STUDY, IS PART OF THE LONG AND CACERES DATASET, AND IS NOT AN RBP.
    labels = [str(x) for x in range(60, 100, 5)]
    num_lc = [51, 50, 49, 43, 38, 26, 16, 8]
    num_new = [288, 134, 83, 46, 24, 16, 9, 1]
        
    plt.scatter(num_lc, num_new)
    plt.plot([52, 52], [-10, 320], linestyle='--', color='#d62728')
    for i in range(len(labels)):
        plt.text(num_lc[i]-5, num_new[i], labels[i]+'%', fontname='Arial', fontsize=13)
        
    plt.xticks(fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel('# of Known SR/SR-related\nProteins Detected', fontname='Arial', fontsize=14)
    plt.ylabel('# of New SR/SR-related\nProtein Candidates', fontname='Arial', fontsize=14)
    
    plt.xlim(0, 57)
    plt.ylim(-2, 310)
    
    fig = plt.gcf()
    fig.set_size_inches(6,3.5)
    plt.savefig('Fig S2B - VariedCombinedSRcomposition_NewSRprots_vs_LCprotsDetected.tiff', bbox_inches='tight', dpi=600)
    plt.close()

if __name__ == '__main__':
    main()