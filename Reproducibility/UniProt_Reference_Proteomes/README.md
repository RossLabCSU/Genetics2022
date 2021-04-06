# Reproducing identification and Pfam analyses of RS domains among UniProt reference proteomes

### Instructions
1. Download all files in the UniProt_Reference_Proteomes directory.
2. Download all archaeal, bacterial, eukaryotic, and viral proteomes from https://figshare.com/articles/dataset/Archaea_ProteinSequences/12937637, (https://figshare.com/articles/dataset/Bacteria_ProteinSequences_Part1/12939812, https://figshare.com/articles/dataset/Bacteria_ProteinSequences_Part2/12939980), (https://figshare.com/articles/dataset/Eukaryota_ProteinSequences_Part1/12937703, https://figshare.com/articles/dataset/Eukaryota_ProteinSequences_Part2/12939479), and https://figshare.com/articles/dataset/Viruses_ProteinSequences/12942362 respectively. NOTE: the bacterial and eukaryotic sets were split into two downloads due to file size constraints.
3. Extract the proteome files *__in separate folders corresponding to each domain of life (i.e. one folder each for Archaea, Bacteria, Eukaryota, and Viruses)__*. In each of these folders, copy the LCD-Composer_MultiProteomeVersion.py script and the appropriate script to make each batch file (ending in "BatchFile.py"). By default, the LCD-Composer_MultiProteomeVersion.py script will scan all sub-folders for FASTA files, so each of the folders should only contain the above-mentioned files with no other subfolders or FASTA files.
4. Navigate to appropriate folder(s) via command line.
5. Run the following commands (NOTE: the commands corresponding to different domains of life can be run concurrently to decrease overall computation time. Additionally, batch files can be divided into smaller batch files and run concurrently to further decrease computation time.):

In the Archaea folder:
```    
python Run_GOanalyses.py
```

```    
python get_GO_Annotations_Q8N5F7.py
```

These two commands generate data appearing in Fig 1D, Fig S3B, Fig S3D, Table S7, and Table S8.
