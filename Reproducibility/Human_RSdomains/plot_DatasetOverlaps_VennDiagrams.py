
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

def main():

    sr_prots = get_SRprots()
    lc_prots = get_protset('All_SRprots_Long_and_Caceres_2009')
    rbps = get_protset('Complete_list_of_RBPs')
    new_SR_RBPs = get_new_SR_RBPs(sr_prots, rbps, lc_prots)

    datasets = ['MixedChargeProts_PhosphorylationConsidered', 'Calarco_et-al_TableS1_HumanHomologs']
    dataset_labels = ['Greig', 'Calarco']
    
    for i in range(len(datasets)):
        dataset = datasets[i]
        dataset_label = dataset_labels[i]
        ds_prots = get_protset(dataset)
                
        sr_only, lc_only, lc_and_sr, ds_only, sr_and_ds, lc_and_ds, lc_and_sr_and_ds = get_dataset_overlaps(sr_prots, rbps, lc_prots, ds_prots, dataset_label)

        if 'MixedCharge' in dataset:
            venn = venn3(subsets = (sr_only, lc_only, lc_and_sr, ds_only, sr_and_ds, lc_and_ds, lc_and_sr_and_ds), set_labels=('My SR Proteins', 'Known SR Proteins', 'Greig et al. Proteins'), alpha=0.5)
        else:
            venn = venn3(subsets = (sr_only, lc_only, lc_and_sr, ds_only, sr_and_ds, lc_and_ds, lc_and_sr_and_ds), set_labels=('My SR Proteins', 'Known SR Proteins', 'Calarco et al. Proteins'), alpha=0.5)
        
        for text in venn.set_labels:
            text.set_fontsize(16)

        fig = plt.gcf()
        fig.set_size_inches(5, 5)
        
        if 'MixedCharge' in dataset:
            plt.savefig('Fig S3A - MixedChargeProt_Comparison_VennDiagram.tiff', bbox_inches='tight', dpi=600)
        else:
            plt.savefig('Fig S3C - MouseHomolog_Comparison_VennDiagram.tiff', bbox_inches='tight', dpi=600)
        plt.close()
        
        
def get_dataset_overlaps(sr_prots, rbps, lc_prots, ds_prots, dataset_label):
    sr_only = [x for x in sr_prots if x not in rbps and x not in lc_prots and x not in ds_prots]
    lc_only = [x for x in lc_prots if x not in rbps and x not in sr_prots and x not in ds_prots]
    ds_only = [x for x in ds_prots if x not in rbps and x not in lc_prots and x not in sr_prots]
    lc_and_ds_notRBP_notSR = [x for x in ds_prots if x in lc_prots and x not in sr_prots]
    sr_and_ds_and_rbp_notLC = [x for x in sr_prots if x in rbps and x not in lc_prots and x in ds_prots]
    sr_and_ds_and_lc_notRBP = [x for x in sr_prots if x in ds_prots and x in lc_prots and x not in rbps]
    sr_and_ds_notRBP_notLC = [x for x in sr_prots if x not in rbps and x not in lc_prots and x in ds_prots]
    sr_and_lc_and_rbp_notDS = [x for x in sr_prots if x in rbps and x in lc_prots and x not in ds_prots]
    sr_and_lc_notRBP_notDS = [x for x in sr_prots if x not in rbps and x in lc_prots and x not in ds_prots]
    sr_and_rbp_notLC_notDS = [x for x in sr_prots if x in rbps and x not in lc_prots and x not in ds_prots]
    lc_and_rbp_notSR_notDS = [x for x in lc_prots if x in rbps and x not in sr_prots and x not in ds_prots]
    ds_and_rbp_notSR_notLC = [x for x in ds_prots if x in rbps and x not in sr_prots and x not in lc_prots]
    sr_andRBP_andLC_andDS = [x for x in sr_prots if x in rbps and x in lc_prots and x in ds_prots]
    
    print('\n\n================================================================================')
    print(dataset_label)
    print('____________')
    print('KEY: \n\t"SR"=proteins identified in this study\n\t"LC"=proteins defined in Long and Caceres (2009)\n\t"RBP"=Known RNA-binding proteins\n\t"Greig"=third dataset of proteins representing mixed-charge proteins, defined in Greig et al. 2020\n\t"Calarco"=third dataset of proteins representing human homologs of mouse SR/SR-related proteins, defined in Calarco et al. 2009)\n')
    print('SR only:', len(sr_only))
    print('LC only:', len(lc_only))
    print(dataset_label + ' only:', len(ds_only))
    print('LC+'+dataset_label+':', len(lc_and_ds_notRBP_notSR))
    print('SR+'+dataset_label+'+RBP:', len(sr_and_ds_and_rbp_notLC))
    print('SR+'+dataset_label+':', len(sr_and_ds_notRBP_notLC))
    print('SR+LC+RBP:', len(sr_and_lc_and_rbp_notDS))
    print('SR+LC+'+dataset_label+':', len(sr_and_ds_and_lc_notRBP))
    print('SR+LC:', len(sr_and_lc_notRBP_notDS))
    print('SR+RBP:', len(sr_and_rbp_notLC_notDS))
    print('LC+RBP:', len(lc_and_rbp_notSR_notDS))
    print(dataset_label+'+RBP:', len(ds_and_rbp_notSR_notLC))
    print('SR+RBP+LC+'+dataset_label+':', len(sr_andRBP_andLC_andDS))
    print('\n\n================================================================================')
    
    # FOR VENN DIAGRAM============================================================================
    sr_only = [x for x in sr_prots if x not in lc_prots and x not in ds_prots]
    lc_only = [x for x in lc_prots if x not in sr_prots and x not in ds_prots]
    ds_only = [x for x in ds_prots if x not in lc_prots and x not in sr_prots]
    lc_and_sr = [x for x in sr_prots if x in lc_prots and x not in ds_prots]
    lc_and_ds = [x for x in lc_prots if x in ds_prots and x not in sr_prots]
    sr_and_ds = [x for x in sr_prots if x in ds_prots and x not in lc_prots]
    sr_and_ds = [x for x in sr_prots if x in ds_prots and x not in lc_prots]
    lc_and_sr_and_ds = [x for x in sr_prots if x in ds_prots and x in lc_prots]

    return [len(x) for x in [sr_only, lc_only, lc_and_sr, ds_only, sr_and_ds, lc_and_ds, lc_and_sr_and_ds]]


def get_protset(dataset):

    h = open(dataset + '.txt')
    prots = set()
    for line in h:
        prots.add(line.rstrip())
    h.close()
    
    return list(prots)

        
def get_SRprots():
    
    h = open('TableS1_Human_RSdomains.tsv')
    header = h.readline()
    
    prots = []
    for line in h:
        items = line.rstrip().split('\t')
        prots.append(items[0])
        
    h.close()
    
    return prots
    
    
def get_new_SR_RBPs(sr_prots, rbps, lc_prots):

    prots = []
    for prot in sr_prots:
        if prot in rbps and prot not in lc_prots:
            prots.append(prot)
            
    return prots


if __name__ == '__main__':
    main()