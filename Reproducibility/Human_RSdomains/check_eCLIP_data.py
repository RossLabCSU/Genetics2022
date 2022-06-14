
import statistics
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

def main():

    encode_pairs = get_ENCODEprots()
    file = 'catrapid_human_basic_normalised.txt'
    num_lines = get_filelength(file)
    h = open(file)
    output = open('catrapid_human_basic_normalised_ENCODEprots_Only_BlueCurve_HighConfidenceOnly.txt', 'w')
    header = h.readline()
    output.write(header)
    count = 0
    
    scores = []
    for line in tqdm(h, total=num_lines):
        prot, transcript, score = line.rstrip().split('\t')
        if (prot, transcript) not in encode_pairs:
            continue
            
        score = float(score)
        scores.append(score)
        
        output.write(line)
        count += 1
        
    print('# of protein-RNA interactions: ', count)
    
    ave = statistics.mean(scores)
    stddev = statistics.stdev(scores)
    print('Mean:', ave, 'Standard deviation:', stddev, 'Z-score=1 @ ' + str(ave+stddev))
    histogram(scores)
    
    
def histogram(scores):

    sns.distplot(scores)
    plt.savefig('eCLIP high confidence interaction catRAPID scores.jpg', bbox_inches='tight', dpi=600)
    plt.close()
    
    
def get_filelength(file):
    
    h = open(file)
    count = 0
    for line in h:
        count += 1
        
    h.close()
    
    return count
    

def get_ENCODEprots():

    h = open('eclip.txt')
    h.readline()
    
    pairs = set()
    for line in h:
        items = line.rstrip().split('\t')
        hc = items[-1]
        species = items[5]
        if species != 'human':
            print(line)
        if int(hc) < 1:
            continue
        prot = items[1]
        cell_lines = items[6]
        transcript = items[2]
        pairs.add( (prot, transcript) )
        
    h.close()
    
    print('Length of ENCODE prots:', len(pairs))
    
    return pairs

if __name__ == '__main__':
    main()