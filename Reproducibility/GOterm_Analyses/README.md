# Reproducing human SR/SR-related protein GO term analyses

### Instructions
1. Download all files in the GOterm_Analyses directory.
2. Download the "HumanReferenceProteom.zip and the "All_SRprots_Long_and_Caceres_2009.txt" file from the Human_RSdomains directory and place them in the same location as the files from Step 1.
3. Download the "LCD-Composer_MultiProteome_MaxCompThreshold.py" script, the "TaxonID_to_Organism.zip" file, and the "TaxonID_to_CommonOrganismName.dat" file from the UniProt_Reference_Proteomes directory and place them in the same location as the files from Steps 1 and 2.
4. Extract all compressed files in the same location as the downloaded files.
5. Navigate to appropriate folder via command line.
6. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):

```
python LCD-Composer_MultiProteome_MaxCompThreshold.py Human_R-richOnly_Results -a R -c 35 -x R_S -m 40_5
```

```
python LCD-Composer_MultiProteome_MaxCompThreshold.py Human_S-richOnly_Results -a S -c 35 -x S_R -m 40_5
```

```
python gather_S-richOnly_and_R-richOnly_prots.py
```

```    
python Run_GOanalyses.py
```

```    
python get_GO_Annotations_Q8N5F7.py
```

This series of commands generates data appearing in Fig 1D, Fig S3B, Fig S3D, Tables S2-S4, Table S11, and Table S12.
