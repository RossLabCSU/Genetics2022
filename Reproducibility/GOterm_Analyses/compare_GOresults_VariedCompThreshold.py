
from scipy import stats
import math

def main():

    files = ['Fig1D_NewRBPs_LCprotsExcluded.tsv', 'TableSY_NewSR-related_RBPs_65percCompThreshold_GOresults.tsv', 'TableSX_NewSR-related_RBPs_60percCompThreshold_GOresults.tsv']
    go_set = get_relevant_goterms(files)
    
    threshold70_df = get_enriched_data('Fig1D_NewRBPs_LCprotsExcluded.tsv', go_set)
    threshold65_df = get_enriched_data('TableSY_NewSR-related_RBPs_65percCompThreshold_GOresults.tsv', go_set)
    threshold60_df = get_enriched_data('TableSX_NewSR-related_RBPs_60percCompThreshold_GOresults.tsv', go_set)
    cats = ['Ratio in Study', 'Ratio in Population', 'lnOR', 'Sidak p-value', 'Associated Proteins']
    output = open('Table S6 - GOtermComparison_VariedCompThreshold.tsv', 'w')
    output.write('GO term\t' + '\t'.join( ['70% S+R, ' + cat + '\t' + '65% S+R, ' + cat + '\t' + '60% S+R, ' + cat for cat in cats] ) + '\n')

    for goterm in go_set:
        output.write(goterm)
        for cat in cats:
            for df in [threshold70_df, threshold65_df, threshold60_df]:
                if goterm in df:
                    output.write( '\t' + str(df[goterm][cat]) )
                else:
                    output.write( '\tN/A' )
                    
        output.write('\n')

    output.close()
            
      
def get_enriched_data(file, go_set):

    df = {}
    h = open(file)
    header = h.readline()
    
    for line in h:
        items = line.rstrip().split('\t')
        go_term = items[3]
        if go_term not in go_set:
            continue
            
        sidak_pval = float(items[10])

        go_term, ratio_in_study, ratio_in_pop = items[3:6]
        proteins = items[12]
        sr_hits, sr_totalprots = ratio_in_study.split('/')
        proteome_hits, proteome_totalprots = ratio_in_pop.split('/')
        sr_hits, sr_totalprots, proteome_hits, proteome_totalprots = [int(x) for x in [sr_hits, sr_totalprots, proteome_hits, proteome_totalprots]]
        sr_totalprots = sr_totalprots - sr_hits
        nonsr_hits = proteome_hits - sr_hits
        nonsr_totalprots = proteome_totalprots - sr_totalprots - proteome_hits
        contingency = [[sr_hits, sr_totalprots], [nonsr_hits, nonsr_totalprots]]
        oddsratio, pval = stats.fisher_exact(contingency, alternative='two-sided')
        lnOR = math.log(oddsratio)
        
        df[go_term] = {'Ratio in Study':ratio_in_study,
                    'Ratio in Population':ratio_in_pop,
                    'lnOR':lnOR,
                    'Sidak p-value':sidak_pval,
                    'Associated Proteins':proteins}
                    
    h.close()
    
    return df
    
    
def get_relevant_goterms(files):
    
    go_set = set()
    for file in files:
        h = open(file)
        header = h.readline()
        
        for line in h:
            items = line.rstrip().split('\t')
            sidak_pval = float(items[10])
            if sidak_pval >= 0.05 or items[2] == 'p':
                continue
                
            go_term = items[3]
            go_set.add(go_term)
        h.close()
        
    return go_set
        
    
if __name__ == '__main__':
    main()