
from scipy import stats
import mpmath
mpmath.mp.dps = 300
import math
import pickle
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']

def main():

    sr_pfam_df = pickle.load(open('SR_proteins_All_Pfams_dict_SHARED_PROTS_EXCLUDED.dat', 'rb'))
    sr_top10_pfams = pickle.load(open('SR_proteins_top10_Pfams_dict.dat', 'rb'))
    residues = ['R', 'S']
    
    output = open('TableS10_Pfam_Enrichment_SRprots.tsv', 'w')
    output.write('\t'.join( ['Domain of Life', 'Pfam Annotation', 'Comparison Group', '# of Occurrences in SR/SR-related Proteins', '# of All Other Pfams in SR/SR-related Proteins', '# of Occurrences in Comparsion Group', '# of All Other Pfams in Comparison Group', 'Odds Ratio', 'lnOR', '95% CI Lower Bound for lnOR Estimate', '95% CI Upper Bound for lnOR Estimate', 'P-value', 'Holm-Sidak Corrected P-value'] ))
    output.write('\t\tNOTE: positive lnOR values indicate enrichment of the associated Pfam annotation in the SR/SR-related protein set relative to the "Comparison Group", whereas negative lnOR values indicate enrichment of the associated Pfam annotation in the "Comparison Group" relative to the "Protein Set of Origin"\n')
    
    for i in range(len(residues)):
        aa = residues[i]
        if aa == 'R':
            limited_aa = 'S'
        else:
            limited_aa = 'R'
            
        other_pfam_df = pickle.load(open(aa+'richOnly_proteins_All_Pfams_dict_SHARED_PROTS_EXCLUDED.dat', 'rb'))

        for domain in domains:
            stored_output_strings = []
            pvals = []
            for pfam in set(sr_top10_pfams[domain]):
                contingency, hits_sr, nonhits_sr, hits_other, nonhits_other = make_contingency(pfam, domain, sr_pfam_df, other_pfam_df)
                oddsratio, lnOR, pval, upper_CI, lower_CI = calc_logodds(contingency, hits_sr, nonhits_sr, hits_other, nonhits_other)
                if str(pval) == '0.0':
                    output_str = '\t'.join( [str(x) for x in [domain, pfam, aa+'rich Proteins', hits_sr, nonhits_sr, hits_other, nonhits_other, oddsratio, lnOR, lower_CI, upper_CI, '<1.0E-300']] )
                else:
                    output_str = '\t'.join( [str(x) for x in [domain, pfam, aa+'rich Proteins', hits_sr, nonhits_sr, hits_other, nonhits_other, oddsratio, lnOR, lower_CI, upper_CI, pval]] )
                pvals.append(pval)
                stored_output_strings.append(output_str)

            corrected_pvals = sidak_correction(pvals)
            for i in range(len(stored_output_strings)):
                corrected_pval = corrected_pvals[i]
                if corrected_pval == '0.0':
                    output.write(stored_output_strings[i] + '\t<1.0E-300\n')
                else:
                    output.write(stored_output_strings[i] + '\t' + str(corrected_pvals[i]) + '\n')
                
    output.close()

    
def sidak_correction(pvals):
    
    positions = [i for i in range(len(pvals))]
    sorted_pvals, pval_positions = zip(*sorted(zip(pvals, positions)))
    m = len(pvals)

    mpmath.mp.dps = 300

    corrected_pvals = [1-(1-mpmath.mp.mpf(sorted_pvals[len(pvals)-i]))**(i) for i in range(m, 0, -1)]
    
    final_pvals = [corrected_pvals[0]]
    for i in range(1, len(corrected_pvals)):
        pval = max(final_pvals[-1], corrected_pvals[i])
        final_pvals.append(pval)
        
    corrected_pvals = [corrected_pvals[0]] + [max(corrected_pvals[i-1], corrected_pvals[i]) for i in range(1, len(corrected_pvals))]
    corrected_pvals = [mpmath.nstr(x, 16) for x in corrected_pvals]
    final_pvals = [mpmath.nstr(x, 16) for x in final_pvals]

    pval_positions, final_pvals = zip(*sorted(zip(pval_positions, final_pvals)))

    return final_pvals
    
            
def make_contingency(pfam, domain, sr_pfam_df, other_pfam_df):
    
    hits_other = other_pfam_df[domain].count(pfam)
    nonhits_other = len(other_pfam_df[domain]) - hits_other
    hits_sr = sr_pfam_df[domain].count(pfam)
    nonhits_sr = len(sr_pfam_df[domain]) - hits_sr

    contingency = [[hits_sr, nonhits_sr], [hits_other, nonhits_other]]  # THIS VERSION TESTS FOR ENRICHMENT AMONG SR PROTEIN SET

    return contingency, hits_sr, nonhits_sr, hits_other, nonhits_other
    
    
def calc_logodds(contingency, hits_sr, nonhits_sr, hits_other, nonhits_other):
    
    oddsratio, pval = stats.fisher_exact(contingency, alternative='two-sided')
    
    if math.isinf(oddsratio):
        oddsratio = 'N/A'
        lnOR = 'N/A'
        upper_CI = 'N/A'
        lower_CI = 'N/A'
    else:    
    
        lnOR = math.log(oddsratio)
            
        upper_CI = lnOR + 1.96 * math.sqrt(1/hits_other + 1/nonhits_other + 1/hits_sr + 1/nonhits_sr)   #95% UPPER CI FOR THE lnOR
        lower_CI = lnOR - 1.96 * math.sqrt(1/hits_other + 1/nonhits_other + 1/hits_sr + 1/nonhits_sr)   #95% LOWER CI FOR THE lnOR
    
    return oddsratio, lnOR, pval, upper_CI, lower_CI


if __name__ == '__main__':
    main()