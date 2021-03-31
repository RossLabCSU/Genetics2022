

def main():

    output = open('Hsapiens_RNAbinding_Proteins_GOid0003723.txt', 'w')
        
    rbps = set()
    
    h = open('Hsapiens.gaf')
    for line in h:
        if line.startswith('!'):
            continue
        items = line.rstrip().split('\t')
        prot = items[1]
        go_term = items[4]
        if go_term == 'GO:0003723':
            rbps.add(prot)
            
    for prot in rbps:
        output.write(prot + '\n')
        
    h.close()
    output.close()

if __name__ == '__main__':
    main()