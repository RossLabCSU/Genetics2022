
import pickle
import os

def main():

    df = {}
    domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
    index = 0
    for folder in ['./' + domain for domain in domains]:
        domain = domains[index]
        df[domain] = []
        for (dirname, dir, files) in os.walk(folder):
            for file in files:
                if file.endswith('.fasta') and 'UP' in file:
                    fname, ext = file.split('.')
                    if 'additional' in fname:   # SKIPS THE "additional" FILES THAT EXIST FOR SOME PROTEOMES AND CONTAIN ADDITIONAL PROTEIN ISOFORMS. THE PURPOSE HERE IS TO ONLY PLOT THE FREQUENCY DATA FOR THE ONE-PROTEIN-PER-GENE VERSION OF EACH PROTEOME.
                        continue
                    df[domain].append(fname)
        index += 1

    pf = open('DomainsOfLife_AllOrganism_Filenames.dat', 'wb')
    pickle.dump(df, pf)
    pf.close()
    
    
if __name__ == '__main__':
    main()