
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

def main():

    combined_thresholds = [60, 65, 70, 75, 80, 85, 90, 95]
    lc_prots = get_LongCaceres_prots()
    rbps = get_rbps()
        
    for threshold in combined_thresholds:
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
                if sr_prots[prot] == '1':
                    lc_and_sr_and_rbp += 1
            
            if prot in rbps and prot not in sr_prots:
                lc_and_rbp += 1
                miss_prots.append(prot)
                
            if prot not in rbps and prot not in sr_prots:
                lc_only += 1
                
        sr_only = 0
        sr_and_rbp = 0
        for prot in sr_prots:
            if prot not in lc_prots and sr_prots[prot] == '0':
                sr_only += 1
            if prot in rbps and prot not in lc_prots:
                sr_and_rbp += 1
        
        # HARD-CODED 0 IN THE SUBSETS LIST REPRESENTS THE RBPs GROUP ALONE. THEREFORE, ONLY RBPS FOUND FOUND IN EITHER OUR NEW SR PROTEIN SET OR THE KNOWN SR PROTEIN SET ARE CONSIDERED.
        venn = venn3(subsets = (sr_only, lc_only, lc_and_sr, 0, sr_and_rbp, lc_and_rbp, lc_and_sr_and_rbp), set_labels=('', '', ''), alpha=0.5)

        for text in venn.set_labels:
            text.set_fontsize(16)

        fig = plt.gcf()
        fig.set_size_inches(5, 5)
        
        if threshold == 70:
            plt.savefig('Fig 1C - Human_RSdomainProteins_Combined-SR-Threshold_' + str(threshold) + '.tiff', bbox_inches='tight', dpi=600)
        else:
            plt.savefig('Fig S2 - Human_RSdomainProteins_Combined-SR-Threshold_' + str(threshold) + '.tiff', bbox_inches='tight', dpi=600)
        plt.close()
        
    
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