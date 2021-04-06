# Reproducing identification and Pfam analyses of RS domains among UniProt reference proteomes

### Instructions
1. Download all files in the UniProt_Reference_Proteomes directory.
2. Download all archaeal, bacterial, eukaryotic, and viral proteomes from https://figshare.com/articles/dataset/Archaea_ProteinSequences/12937637, (https://figshare.com/articles/dataset/Bacteria_ProteinSequences_Part1/12939812, https://figshare.com/articles/dataset/Bacteria_ProteinSequences_Part2/12939980), (https://figshare.com/articles/dataset/Eukaryota_ProteinSequences_Part1/12937703, https://figshare.com/articles/dataset/Eukaryota_ProteinSequences_Part2/12939479), and https://figshare.com/articles/dataset/Viruses_ProteinSequences/12942362 respectively. NOTE: the bacterial and eukaryotic sets were split into two downloads due to file size constraints.
3. Extract all compressed files in the same location as the downloaded files and scripts.
4. Navigate to appropriate folder via command line.
5. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):

```    
python Run_GOanalyses.py
```

```    
python get_GO_Annotations_Q8N5F7.py
```

These two commands generate data appearing in Fig 1D, Fig S3B, Fig S3D, Table S7, and Table S8.
