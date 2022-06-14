
import tqdm

def main():

    file = 'catrapid_human_basic_normalised.txt'
    num_lines = get_filelength(file)
    h = open(file)
    h.readline()

    output = open('catrapid_human_basic_normalised_AllProteins_SCORE_ATLEAST_21p297.txt', 'w')
    output.write('\t'.join(['uniprot_accession', 'ensembl_transcript_id', 'catrapid_score_max_normalised']) + '\n')
    
    rbps = set()

    for line in tqdm.tqdm(h, total=num_lines):
        prot, transcript, score = line.rstrip().split('\t')
        score = float(score)

        if score >= 21.297:
            output.write(line)
        
    h.close()
    output.close()
    
    
def get_filelength(file):
    
    h = open(file)
    count = 0
    for line in h:
        count += 1
        
    h.close()
    
    return count


if __name__ == '__main__':
    main()