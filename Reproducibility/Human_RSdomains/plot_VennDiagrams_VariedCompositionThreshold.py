
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

def main():

    combined_thresholds = [60, 65, 70, 75, 80, 85, 90, 95]
    lc_prots = get_LongCaceres_prots()
    rbps = get_rbps()
    output = open('Table S5 - Human SR and SR-related Proteins by S-R Composition Threshold.tsv', 'w')
    categories = ['lc_and_sr_and_rbp', 'lc_and_sr', 'sr_and_rbp', 'sr_only', 'lc_and_rbp', 'lc_only']
    cat_labels = ['This study + Long/Caceres + Known RBP', 'This study + Long/Caceres', 'This study + Known RBP', 'This study only', 'Long/Caceres + Known RBP', 'Long/Caceres only']
    header_labels = ['# of Identified Proteins, ' + label + '\t' + 'UniProt IDs, ' + label for label in cat_labels]
    output.write('\t'.join( ['Combined S+R Composition Threshold', 'Total # of Long/Caceres Proteins Detected', 'Total # of New SR-related Protein Candidates from This Study'] + header_labels ) + '\n')
    
    for threshold in combined_thresholds:
        
        threshold_df = {cat:set() for cat in categories}
        sr_prots = get_SRprots(threshold, rbps)
        
        lc_and_sr = 0
        lc_and_sr_and_rbp = 0
        lc_and_rbp = 0
        lc_only = 0
        miss_prots = []
        for prot in lc_prots:
            if prot in sr_prots:
                if prot not in rbps:
                    lc_and_sr += 1
                    threshold_df['lc_and_sr'].add(prot)
                if sr_prots[prot] == '1':
                    lc_and_sr_and_rbp += 1
                    threshold_df['lc_and_sr_and_rbp'].add(prot)
            
            if prot in rbps and prot not in sr_prots:
                lc_and_rbp += 1
                miss_prots.append(prot)
                threshold_df['lc_and_rbp'].add(prot)
                
            if prot not in rbps and prot not in sr_prots:
                lc_only += 1
                threshold_df['lc_only'].add(prot)
                
        sr_only = 0
        sr_and_rbp = 0
        for prot in sr_prots:
            if prot not in lc_prots and sr_prots[prot] == '0':
                sr_only += 1
                threshold_df['sr_only'].add(prot)
            if prot in rbps and prot not in lc_prots:
                sr_and_rbp += 1
                threshold_df['sr_and_rbp'].add(prot)
        
        # HARD-CODED 0 IN THE SUBSETS LIST REPRESENTS THE RBPs GROUP ALONE. THEREFORE, ONLY RBPS FOUND FOUND IN EITHER OUR NEW SR PROTEIN SET OR THE KNOWN SR PROTEIN SET ARE CONSIDERED.
        venn = venn3(subsets = (sr_only, lc_only, lc_and_sr, 0, sr_and_rbp, lc_and_rbp, lc_and_sr_and_rbp), set_labels=('', '', ''), alpha=0.5)

        for text in venn.set_labels:
            text.set_fontsize(16)

        fig = plt.gcf()
        fig.set_size_inches(5, 5)
        
        if threshold == 70:
            plt.savefig('Fig 1C - Human_RSdomainProteins_Combined-SR-Threshold_' + str(threshold) + '.tiff', bbox_inches='tight', dpi=600)
        else:
            plt.savefig('Fig S2A - Human_RSdomainProteins_Combined-SR-Threshold_' + str(threshold) + '.tiff', bbox_inches='tight', dpi=600)
        plt.close()
        
        output.write('\t'.join([str(threshold), str(sum([len(threshold_df['lc_and_sr']), len(threshold_df['lc_and_sr_and_rbp'])])), str(sum([len(threshold_df['sr_only']), len(threshold_df['sr_and_rbp'])]))]))
        for cat in categories:
            output.write('\t' + str(len(threshold_df[cat])))
            output.write('\t' + ', '.join(threshold_df[cat]))
        output.write('\n')
        
    output.close()
        
        
    
def get_SRprots(threshold, rbps):

    df = {}
    for s_comp in range(20, 105, 5):
        for r_comp in range(20, 105, 5):
            if s_comp + r_comp > 100:
                continue

            h = open('Human_S-R_' + str(s_comp) + '-' + str(r_comp) + '_RESULTS.tsv')
            for i in range(7):
                h.readline()
            header = h.readline()
            for line in h:
                items = line.rstrip().split('\t')
                junk, uniprot, *junk = items[0].split('|')

                if s_comp + r_comp >= threshold:
                    if uniprot in rbps:
                        is_rbp = '1'
                    else:
                        is_rbp = '0'
                    df[uniprot] = is_rbp

            h.close()
    
    return df
    
        
def get_LongCaceres_prots():

    h = open('All_SRprots_Long_and_Caceres_2009.txt')
    prots = []
    for line in h:
        prots.append(line.rstrip())
    h.close()
        
    return prots
    
    
def get_rbps():
    
    h = open('Complete_list_of_RBPs.txt')
    rbps = []
    for line in h:
        prot = line.rstrip()
        rbps.append(prot)
    h.close()
    
    return rbps


if __name__ == '__main__':
    main()
