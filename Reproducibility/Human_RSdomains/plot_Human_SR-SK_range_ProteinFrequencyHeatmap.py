
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def main():

    for dataset in ['S-R', 'S-K']:
    
        pos_aa = dataset[-1]
        matrix = []
        
        for s_comp in range(20, 105, 5):
            col = []
            for pos_comp in range(20, 105, 5):
                if s_comp + pos_comp > 100:
                    col.append(0.01)    # DUMMY VALUE USED FOR MASKING CELLS
                    continue

                total_lcds = get_total_lcds(s_comp, pos_comp, dataset)
                if total_lcds == 0:
                    col.append(0.1) # DUMMY VALUE FOR COLOR SCALE ONLY: LATER CONVERTED TO 0
                    continue

                col.append(total_lcds)
                
            matrix.append(col)

        plot_heatmap(matrix, dataset, pos_aa)
    
    
def plot_heatmap(matrix, dataset, pos_aa):

    matrix = np.matrix.transpose(np.array(matrix))
    mask = np.where(matrix<0.05, True, False)     # MASK CELLS WITH S+R COMPOSITION THRESHOLD EXCEEDING 100
    
    log_mat = np.log10(matrix)  # LOG TRANSFORM MATRIX FOR HEATMAP CELL COLORS
    log_mat = np.where(log_mat<0, 0, log_mat)   # CONVERT DUMMY VALUES TO 0
    
    colors = sns.color_palette("coolwarm", 1000)
    
    sns.heatmap(log_mat, cmap=colors, fmt='.0f', annot=matrix, mask=mask)
    plt.xticks([x+0.5 for x in range(17)], labels=[str(x) for x in range(20, 105, 5)])
    plt.yticks([x+0.5 for x in range(17)], labels=[str(x) for x in range(20, 105, 5)], rotation=0)
    plt.xlabel('S Composition Threshold')
    plt.ylabel(pos_aa + ' Composition Threshold')
    plt.ylim(0, 13)
    plt.xlim(0, 13)
    fig = plt.gcf()
    fig.set_size_inches(8,5)
    if dataset == 'S-R':
        plt.savefig('Fig 1A - Human Pairwise ' + dataset + 'comp Range LCD Frequencies.tiff', bbox_inches='tight', dpi=600)
    else:
        plt.savefig('Fig S1A - Human Pairwise ' + dataset + 'comp Range LCD Frequencies.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    
    
def get_total_lcds(s_comp, pos_comp, dataset):

    h = open('Human_' + dataset + '_' + str(s_comp) + '-' + str(pos_comp) + '_RESULTS.tsv')
    for i in range(7):
        h.readline()
    header = h.readline()
    total_lcds = 0
    for line in h:
        total_lcds += 1
    h.close()

    return total_lcds
    
    
if __name__ == '__main__':
    main()
    
