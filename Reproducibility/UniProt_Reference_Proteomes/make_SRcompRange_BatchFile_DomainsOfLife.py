
def main(args):

    domain = args.domain
    output = open('RUN_LCD-Composer_' + domain + '_SRcompRANGE_Batch.bat', 'w')
    
    for s_comp in range(20, 105, 5):
        for r_comp in range(20, 105, 5):
            if s_comp + r_comp > 100:
                continue
            
            output.write('python LCD-Composer_MultiProteome_MaxCompThreshold.py ' + domain + '_S-R_' + str(s_comp) + '-' + str(r_comp) + '_RESULTS -a S_R -c ' + str(s_comp) + '_' + str(r_comp) + '\n')

    output.close() 
    
    
def get_args(arguments):

    parser = argparse.ArgumentParser(description='Identification of low-complexity domains on the basis of amino acid composition and linear dispersion', prog='LCD-Composer')
    parser.add_argument('domain', help="""Domain of life. Should be "Archaea", "Bacteria", "Eukaryota", or "Viruses".""")
    args = parser.parse_args(arguments)
    
    return args
    

if __name__ == '__main__':
    import sys, argparse
    args = get_args(sys.argv[1:])
    main(args)