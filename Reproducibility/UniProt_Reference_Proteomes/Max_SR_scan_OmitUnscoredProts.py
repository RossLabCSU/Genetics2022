
import os
from Bio import SeqIO
from tqdm import tqdm
import datetime
import sys
import argparse

def main(args):

    domain = args.domain_of_life
    coronavirus_proteomes = None
    if domain == 'Coronaviruses':
        coronavirus_proteomes = get_coronavirus_proteomes()
    perc_threshold = args.percent_threshold
    win_sizes = [20, 30, 40, 50]
    min_counts = [int(x*perc_threshold/100) for x in win_sizes]
    
    output = open('Max_SR_scan_' + domain + '_Results_' + str(perc_threshold) + 'percMin_UnscoredProtsOmitted.tsv', 'w')
    output.write('\t'.join(['File of Origin', 'UniProt ID'] + ['Maximum S+R Percentage in ' + str(win_size) + 'aa Window' for win_size in (20, 30, 40, 50)]) + '\n')
    print(domain, str(datetime.datetime.now()))
    
    for dirname, dirs, files in os.walk('.'):
        for file in tqdm(files):
            if not file.startswith('UP0') or not file.endswith('.fasta'):
                continue
                
            if domain == 'Coronaviruses' and file not in coronavirus_proteomes:
                continue

            h = open(file)
            for seq_record in SeqIO.parse(h, 'fasta'):
                prot_id = str(seq_record.id)
                junk, uniprot, *junk = prot_id.split('|')
                seq = str(seq_record.seq)
                if seq[-1] == '*':
                    seq = seq[:-1]
                    
                max_percents = []
                    
                for i, win_size in enumerate(win_sizes):
                    if len(seq) < win_size:
                        max_percents.append('-1')
                        continue
                    min_count = min_counts[i]
                    max_sr = 0
                    for i in range(0, len(seq)-win_size-1):
                        window = seq[i:i+win_size]
                        r_count = window.count('R')
                        s_count = window.count('S')
                        if r_count < min_count or s_count < min_count:
                            continue
                        sr_count = r_count + s_count
                        if sr_count > max_sr:
                            max_sr = sr_count
                            
                    if max_sr == 0:
                        max_sr_perc = -1
                    else:
                        max_sr_perc = max_sr / win_size * 100
                    max_percents.append(str(max_sr_perc))
                    
                # OMIT PROTEINS THAT DO NOT MEET MINIMUM %S AND %R
                if max_percents == ['-1']*4:
                    continue
                    
                output.write('\t'.join([file, uniprot] + max_percents) + '\n')
            h.close()
                    
    output.close()
    
    
def get_coronavirus_proteomes():

    h = open('CoronavirusProteomes_List.txt')
    proteomes = []
    for line in h:
        proteomes.append(line.rstrip())
        
    h.close()
    
    return proteomes
    

def parse_arguments(arguments) :
    
    parser = argparse.ArgumentParser()
    parser.add_argument('domain_of_life', help = 'input the category of life (Archaea, Bacteria, Eukaryota, Viruses, or Coronaviruses)')
    parser.add_argument('-p', '--percent_threshold', type=int, default=10, 
                        help="""Minimum percent composition for S and R (individually) required
                                    for a window to be scored.""")
    args = parser.parse_args(arguments)

    return args

    
if __name__ == '__main__' :
    args = parse_arguments(sys.argv[1:])
    main(args)