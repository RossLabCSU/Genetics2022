
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

pfam_df = pickle.load(open('SelectProtSeqs_Pfam_Results.dat', 'rb'))

def main():

    seqs_df = get_seqs()
    iso_df = get_isoforms()
    iso_df = sort_by_seqlen(iso_df, seqs_df)
    
    plotting(iso_df, seqs_df)
    
    
def get_seqs():

    h = open('SelectProtSeqs.fasta')
    df = {}
    
    for seq_record in SeqIO.parse(h, 'fasta'):
        id = str(seq_record.id)
        seq = str(seq_record.seq)
        junk, prot, *junk = id.split('|')
        df[prot] = seq
    
    h.close()
    
    return df
    
    
def plotting(iso_df, seqs_df):

    colors = sns.color_palette('Set2')
    colors += sns.color_palette('pastel')
    colors += sns.color_palette('dark')
    colors += sns.color_palette('colorblind')
    
    lineEnd_offset = 6
    
    for prot in iso_df:
        for i in range(len(iso_df[prot]['Isoforms'])):
            iso = iso_df[prot]['Isoforms'][i]
            try:
                sr_bounds = iso_df[prot]['SR boundaries'][i]
                sr_bounds = sr_bounds.split(',')
            except:
                sr_bounds = ''
            seq = seqs_df[iso]
            
            sr_positions = []
            all_sr_pos = []
            num_sr_pos = 0
            for bound in sr_bounds:
                if bound == '':
                    continue
                start, end = bound[1:-1].split('-')
                start = int(start)-1 + lineEnd_offset
                end = int(end) - lineEnd_offset
                positions = [x for x in range(start, end)]
                sr_positions.append( positions )
                all_sr_pos += positions
            non_sr_positions = [x for x in range(len(seq)) if x not in all_sr_pos]
            
            plt.plot([x for x in range(len(seq))], [i for aa in range(len(seqs_df[iso]))], color='0.5')
                
            if iso in pfam_df:
                pfams = pfam_df[iso]
            else:
                pfams = []

            color_index = 2
            pfam_label_positions = []
            for pfam in pfams:
                start = int(pfam['env']['from']) - 1 + lineEnd_offset
                end = int(pfam['env']['to']) - lineEnd_offset
                name = pfam['name']
                positions = [x for x in range(start, end)]
                
                if 'RRM' in name:
                    plt.plot(positions, [i for x in range(len(positions))], color=colors[1], linewidth=15, solid_capstyle='round')
                else:
                    plt.plot(positions, [i for x in range(len(positions))], color=colors[color_index], linewidth=15, solid_capstyle='round')
                x_pos = min(positions) + ((max(positions) - min(positions)) / 2)
                if prot == 'SCAF4': # SPECIAL FONT SIZING FOR SCAF4 TO PREVENT TEXT OVERLAP
                    plt.text(x_pos, i, name, ha='center', va='center', fontname='Arial', fontsize=8)
                else:
                    plt.text(x_pos, i, name, ha='center', va='center', fontname='Arial')
                pfam_label_positions.append(x_pos)
                color_index += 1
                
            for sr_domain in sr_positions:
                plt.plot(sr_domain, [i for x in range(len(sr_domain))], color=colors[0], linewidth=15, solid_capstyle='round')
                x_pos = min(sr_domain) + ((max(sr_domain) - min(sr_domain)) / 2)
                plt.text(x_pos, i, 'RS', ha='center', va='center', fontname='Arial')
                
            # MANUALLY ADD THE RS DOMAIN OF SRSF9, INDICATED AS "RS*" BECAUSE IT PASSES A 60% COMBINED S+R THRESHOLD BUT NOT THE STANDARD 70% COMBINED S+R THRESHOLD
            if prot == 'SRSF9' and i == 1:
                bound = '(186-211)'
                start, end = bound[1:-1].split('-')
                start = int(start)-1 + lineEnd_offset
                end = int(end) - lineEnd_offset
                sr_domain = [x for x in range(start, end)]
                
                plt.plot(sr_domain, [i for x in range(len(sr_domain))], color=colors[0], linewidth=15, solid_capstyle='round')
                x_pos = min(sr_domain) + ((max(sr_domain) - min(sr_domain)) / 2)

                plt.text(x_pos, i, 'RS*', ha='center', va='center', fontname='Arial')
                
                
            plt.text(len(seq)+1, i, str(len(seq)), ha='left', va='center', fontname='Arial')
                
        fig = plt.gcf()
        longest_iso = iso_df[prot]['Isoforms'][-1]
        length = len(seqs_df[longest_iso])
        num_isos = len(iso_df[prot]['Isoforms'])
        
        plt.ylim(-1, num_isos+1)
        plt.axis('off')

        if num_isos == 2:
            fig.set_size_inches(7, 1.4)
        elif num_isos == 3:
            fig.set_size_inches(7, 1.8)
        else:
            fig.set_size_inches(7, num_isos/2)
            
        # CAN FAIL FOR SOME PROTEINS IF THEY ARE TOO LONG
        try:
            plt.savefig('FigS6_' + prot + '_IsoformCartoon_FragmentsRemoved.tiff', bbox_inches='tight', dpi=600)
        except:
            fig.set_size_inches(2*(length/200), num_isos)
            plt.savefig('FigS6_' + prot + '_IsoformCartoon.tiff', bbox_inches='tight', dpi=600)
        plt.close()
    
    
def get_isoforms():
    # ALSO COLLECTS UNIPROT IDs FOR ISOFORMS THAT WERE NOT IDENTIFIED AS SR/SR-RELATED PROTEINS
    h = open('TableS5_Human_RSdomains_IsoformComparison.tsv')
    header = h.readline()
    df = {}
    iso_tracker = {}
    for line in h:
        items = line.rstrip().split('\t')
        uniprot = items[0]
        common_name = items[1]
        if 'SRSF' not in common_name and common_name != 'SCAF4' and common_name != 'ZRANB2' and common_name != 'NKAPL':
            continue
            
        bounds = items[13]

        all_iso_uniprots = items[8].split(',')
            
        df[common_name] = df.get(common_name, {'Isoforms':[], 'SR boundaries':[]})
        df[common_name]['Isoforms'].append(uniprot)
        df[common_name]['SR boundaries'].append(bounds)
        iso_tracker[common_name] = all_iso_uniprots
        
    h.close()
    
    for prot in iso_tracker:
        for iso in iso_tracker[prot]:
            if iso not in df[prot]['Isoforms']:
                df[prot]['Isoforms'].append(iso)
                df[prot]['SR boundaries'].append('')
    
    # MANUALLY ADD SRSF9
    df['SRSF9'] = {'Isoforms':[], 'SR boundaries':[]}
    h = open('SelectProtSeqs.fasta')
    for seq_record in SeqIO.parse(h, 'fasta'):
        desc = str(seq_record.description)
        if 'GN=SRSF9' in desc:
            id = str(seq_record.id)
            junk, uniprot, *junk = id.split('|')
            df['SRSF9']['Isoforms'].append(uniprot)
    h.close()
        
    return df
    
    
def sort_by_seqlen(df, seqs_df):
    
    sorted_df = {}
    for prot in df:
        isos = df[prot]['Isoforms']
        bounds = df[prot]['SR boundaries']
        seqlens = [len(seqs_df[acc]) for acc in isos]
        if prot == 'SRSF9':     # SORTING FAILS FOR SRSF9 SINCE THERE ARE NO RS DOMAINS IDENTIFIED IN THIS PROTEIN
            seqlens, isos = zip(*sorted(zip(seqlens, isos)))
        else:
            seqlens, isos, bounds = zip(*sorted(zip(seqlens, isos, bounds)))
            
        sorted_df[prot] = {'Isoforms':isos, 'SR boundaries':bounds}
        
    return sorted_df
        

if __name__ == '__main__':
    main()