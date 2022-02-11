
from Bio import SeqIO

def main():

    output = open('SelectProtSeqs.fasta', 'w')
    
    h = open('UP000005640_9606.fasta')
    for seq_record in SeqIO.parse(h, 'fasta'):
        desc = str(seq_record.description)
        if ('GN=SRSF' in desc or 'GN=SCAF4' in desc or 'GN=ZRANB2' in desc or 'GN=NKAPL' in desc) and '(Fragment)' not in desc and 'readthrough' not in desc:
            seq = str(seq_record.seq)
            output.write('>'+desc+'\n')
            output.write(seq + '\n')
    
    h.close()
    
    h = open('UP000005640_9606_additional.fasta')
    for seq_record in SeqIO.parse(h, 'fasta'):
        desc = str(seq_record.description)
        if ('GN=SRSF' in desc or 'GN=SCAF4' in desc or 'GN=ZRANB2' in desc or 'GN=NKAPL' in desc) and '(Fragment)' not in desc and 'readthrough' not in desc:
            seq = str(seq_record.seq)
            output.write('>'+desc+'\n')
            output.write(seq + '\n')
    
    h.close()
    
    output.close()

if __name__ == '__main__':
    main()