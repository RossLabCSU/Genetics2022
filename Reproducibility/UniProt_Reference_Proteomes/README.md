# Reproducing identification and Pfam analyses of RS domains among UniProt reference proteomes

### Instructions
1. Download all files in the UniProt_Reference_Proteomes directory.
2. Download and extract all of the relevant Pfam annotation files from https://figshare.com/account/articles/19141799.
3. Download all archaeal, bacterial, eukaryotic, and viral proteomes from https://figshare.com/articles/dataset/Archaea_ProteinSequences/12937637, (https://figshare.com/articles/dataset/Bacteria_ProteinSequences_Part1/12939812, https://figshare.com/articles/dataset/Bacteria_ProteinSequences_Part2/12939980), (https://figshare.com/articles/dataset/Eukaryota_ProteinSequences_Part1/12937703, https://figshare.com/articles/dataset/Eukaryota_ProteinSequences_Part2/12939479), and https://figshare.com/articles/dataset/Viruses_ProteinSequences/12942362 respectively. NOTE: the bacterial and eukaryotic sets were split into two downloads due to file size constraints.
4. Extract the proteome files *__in separate folders corresponding to each domain of life. The folders MUST be named "Archaea", "Bacteria", "Eukaryota", and "Viruses", respectively__*.
5. Place a copy of the "LCD-Composer_MultiProteome_MaxCompThreshold.py", "make_SRcompRange_BatchFile_DomainsOfLife.py", "get_SR_prots_above_CombinedSRthreshold_MultiProteome.py", "TaxonID_to_CommonOrganismName.dat", and "TaxonID_to_Organism_ReferenceProteomes.dat" files into each of the folders that you created in Step 4. By default, the LCD-Composer_MultiProteome_MaxCompThreshold.py script will scan all sub-folders for FASTA files, so each of the folders should only contain the above-mentioned files with no other subfolders or FASTA files.
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

7. Once all of the commands listed above have finished running, each folder should contain: 1) a file ending in "SequenceDictionary" (.dat file), 2) a file ending in "SR_proteins_with_Combined_S-R_Above_70" (.tsv file), 3) a file ending in "R-richOnly_RESULTS" (.tsv file), and 4) a file ending in "S-richOnly_RESULTS" (.tsv file). Copy these files from each folder into a single location that also contains the remaining scripts from this Github directory. This location should also contain all of the extracted Pfam annotation files from Step 2 above.
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

These two commands generate Fig 4A and preliminary supplementary tables.
</br></br>
For analyses of Pfam domain annotations among SR/SR-related proteins in the UniProt reference proteomes, run the following commands in-sequence:

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

This series of commands generates Fig 4B-E, as well as the final versions of Tables S8-S12.

For BLAST analyses to identify orthologs of SR-related helicases in archaea and bacteria, you must first:
1. Install BLAST version 2.10.1 on your local machine.
2. Copy all archaeal and bacterial proteomes not ending in "additional.fasta" to the bin folder located in the installed BLAST folder.
3. Copy necessary scripts and text files into the bin folder.
4. Copy TableS8 (Archaea) and TableS9 (Bacteria) generated from the commands above into the bin folder.
5. Run the following commands (IMPORTANT NOTE: in total, these commands will end up generating ~5million files. Opening this folder with a graphical viewer may crash the viewer and possibly your computer. Command-line navigation is recommended):

```
python mask_RSdomains_QuerySequences.py Archaea
```

```
python mask_RSdomains_QuerySequences.py Bacteria
```
