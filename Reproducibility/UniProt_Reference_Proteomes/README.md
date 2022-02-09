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

```
python LCD-Composer_MultiProteome_MaxCompThreshold.py Archaea_R-rich-only_RESULTS -a R -c 35 -x R_S -m 40_5
```

```
python LCD-Composer_MultiProteome_MaxCompThreshold.py Archaea_S-rich-only_RESULTS -a S -c 35 -x S_R -m 40_5
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

```
python LCD-Composer_MultiProteome_MaxCompThreshold.py Bacteria_R-rich-only_RESULTS -a R -c 35 -x R_S -m 40_5
```

```
python LCD-Composer_MultiProteome_MaxCompThreshold.py Bacteria_S-rich-only_RESULTS -a S -c 35 -x S_R -m 40_5
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

```
python LCD-Composer_MultiProteome_MaxCompThreshold.py Eukaryota_R-rich-only_RESULTS -a R -c 35 -x R_S -m 40_5
```

```
python LCD-Composer_MultiProteome_MaxCompThreshold.py Eukaryota_S-rich-only_RESULTS -a S -c 35 -x S_R -m 40_5
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

```
python LCD-Composer_MultiProteome_MaxCompThreshold.py Viruses_R-rich-only_RESULTS -a R -c 35 -x R_S -m 40_5
```

```
python LCD-Composer_MultiProteome_MaxCompThreshold.py Viruses_S-rich-only_RESULTS -a S -c 35 -x S_R -m 40_5
```

7. Once all of the commands listed above have finished running, each folder should contain: 1) a file ending in "SequenceDictionary" (.dat file), 2) a file ending in "SR_proteins_with_Combined_S-R_Above_70" (.tsv file), 3) a file ending in "R-richOnly_RESULTS" (.tsv file), and 4) a file ending in "S-richOnly_RESULTS" (.tsv file). Copy these files from each folder into a single location that also contains the remaining scripts from this Github directory.
8. Navigate to this new location via command line and run the following commands in-sequence:

```    
python Merge_Final_RSdomains.py
```

```    
python get_DomainsOfLife_OrganismList_df.py
```

```    
python plot_NumberOfProts_with_SRdomain_DomainsOfLife.py
```

These two commands generate Fig 4A and Tables S6-S9.
</br></br>
For analyses of Pfam domain annotations among SR/SR-related proteins in the UniProt reference proteomes, run the following commands in-sequence to generate Fig 4B-E, as well as the final versions of Tables S6-S9:

```    
python map_Pfams_to_SRprots.py
```

```    
python map_Pfams_to_Srich_and_Rrich_prots.py
```

```    
python plot_TopPfamAnnotations.py
```

```
python determine_top10_Pfams.py
```

```
python calc_PfamEnrichment_SRprots_vs_Srich_or_Rrich_prots.py
```
