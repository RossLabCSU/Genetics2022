
def main():

    thresholds = [(37,38), (35,37)]
    residues = ['S', 'R']
    other_aas = ['R', 'S']
    
    for i in range(2):
        threshold = thresholds[i]
        residue = residues[i]
        other_aa = other_aas[i]
        
        h = open('Human_'+residue+'-richOnly_Results.tsv')
        prots = set()
        for i in range(8):
            h.readline()
        for line in h:
            items = line.rstrip().split('\t')
            junk, prot, *junk = items[2].split('|')
            if residue == 'R':
                comp = float(items[21])
            else:
                comp = float(items[22])
                
            if comp > threshold[0]+0.00001 and comp < threshold[1]+0.00001:
                prots.add(prot)
                
        output = open('Human_'+residue+'-richOnly_FILTERED_PROTEIN_LIST.txt', 'w')
        output.write('\n'.join(list(prots)))
        output.close()
            
if __name__ == '__main__':
    main()