# Reproducing identification and Pfam analyses of RS domains among UniProt reference proteomes

### Instructions
1. Download all files in the UniProt_Reference_Proteomes directory.
2. Download all archaeal, bacterial, eukaryotic, and viral proteomes from https://figshare.com/articles/dataset/Archaea_ProteinSequences/12937637, (https://figshare.com/articles/dataset/Bacteria_ProteinSequences_Part1/12939812, https://figshare.com/articles/dataset/Bacteria_ProteinSequences_Part2/12939980), (https://figshare.com/articles/dataset/Eukaryota_ProteinSequences_Part1/12937703, https://figshare.com/articles/dataset/Eukaryota_ProteinSequences_Part2/12939479), and https://figshare.com/articles/dataset/Viruses_ProteinSequences/12942362 respectively. NOTE: the bacterial and eukaryotic sets were split into two downloads due to file size constraints.
3. Extract the proteome files *__in separate folders corresponding to each domain of life. The folders MUST be named "Archaea", "Bacteria", "Eukaryota", and "Viruses", respectively__*.
4. Place a copy of the "LCD-Composer_MultiProteome_MaxCompThreshold.py", "make_SRcompRange_BatchFile_DomainsOfLife.py", and "get_SR_prots_above_CombinedSRthreshold_MultiProteome.py" scripts into each of the folders that you created in Step 3. By default, the LCD-Composer_MultiProteome_MaxCompThreshold.py script will scan all sub-folders for FASTA files, so each of the folders should only contain the above-mentioned files with no other subfolders or FASTA files.
5. Download the "TaxonID_to_CommonOrganismName.dat" file and "TaxonID_to_Organism.zip" file from the Human_RSdomain directory, extract the TaxonID_to_Organism file, and place a copy of these files in each of the folders generated in Step 3.
6. Navigate to appropriate folder(s) via command line.
7. Run the following commands (NOTE: the commands corresponding to different domains of life can be run concurrently to decrease overall computation time. Additionally, batch files can be divided into smaller batch files and run concurrently to further decrease computation time.):

In the Archaea folder run the following commands:
```    
python make_SRcompRange_BatchFile_DomainsOfLife.py Archaea
```

```    
.\RUN_LCD-Composer_Archaea_SRcompRANGE_Batch.bat
```

```    
python get_SR_prots_above_CombinedSRthreshold_MultiProteome.py Archaea
```

In the Bacteria folder run the following commands:
```    
python make_SRcompRange_BatchFile_DomainsOfLife.py Bacteria
```

```    
.\RUN_LCD-Composer_Bacteria_SRcompRANGE_Batch.bat
```

```    
python get_SR_prots_above_CombinedSRthreshold_MultiProteome.py Bacteria
```

In the Eukaryota folder run the following commands:
```    
python make_SRcompRange_BatchFile_DomainsOfLife.py Eukaryota
```

```    
.\RUN_LCD-Composer_Eukaryota_SRcompRANGE_Batch.bat
```

```    
python get_SR_prots_above_CombinedSRthreshold_MultiProteome.py Eukaryota
```

In the Viruses folder run the following commands:
```    
python make_SRcompRange_BatchFile_DomainsOfLife.py Viruses
```

```    
.\RUN_LCD-Composer_Viruses_SRcompRANGE_Batch.bat
```

```    
python get_SR_prots_above_CombinedSRthreshold_MultiProteome.py Viruses
```

7. Once all of the commands listed above have finished running, each folder should contain a file ending in "SequenceDictionary" (.fasta file) and a file ending in "SR_proteins_with_Combined_S-R_Above_70" (.tsv file). Copy these files from each folder into a single location that also contains the XXXXXXXXXXXXXXXXXXx scripts.
8. Navigate to this new location via command line and run the following commands in-sequence:

```    
python get_DomainsOfLife_OrganismList_df.py
```

```    
python plot_NumberOfProts_with_SRdomain_DomainsOfLife.py
```

These two commands generate Fig 4A and Tables S3-S6.
</br></br>
For analyses of Pfam domain annotations among SR/SR-related proteins in the UniProt reference proteomes, run the following commands in-sequence:

```    
python fragment_FastaFiles.py
```

```    
python make_Pfam_BatchFile.py
```

```    
.\Archaea_Pfam_BatchFile.bat
```

```    
.\Bacteria_Pfam_BatchFile.bat
```

```    
.\Eukaryota_Pfam_BatchFile.bat
```

```    
.\Viruses_Pfam_BatchFile.bat
```

Note that Pfam server errors sometimes fail to generate results files for a small subset of commands specified in the batch files. In these instances, simply identify the skipped files and re-run the corresponding commands from the batch files.
</br></br>
Run the following commands in-seqeuence to generate Fig 4B-E and Tables S3-S6:

```    
python gather_DomainsOfLife_Pfam_Results.py
```

```    
python map_Pfams_to_SRprots.py
```

```    
python plot_TopPfamAnnotations.py
```
