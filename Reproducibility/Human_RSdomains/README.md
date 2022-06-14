# Reproducing SR and SK searches (and subsequent analyses) for the human proteome

### Instructions
1. Download LCD-Composer.py from https://github.com/RossLabCSU/LCD-Composer, as well as all files in the Human_RSdomains directory.
2. Extract all compressed files in the same location as LCD-Composer.py and other Python scripts. NOTE: the "UP000005640_9606_additional.fasta" file containing all isoform sequences for human proteins has been pre-filtered to exclude sequences that are classified as non-canonical readthrough protein products or as protein fragments.
3. Navigate to appropriate folder via command line.
4. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):

```    
python make_Human_SRsearch_and_SKsearch_BatchFile.py
```

```
.\RUN_LCD-Composer_Human_SRsearch_and_SKsearch_Batch.bat
```

```    
python get_Human_RNAbinding_Prots_from_GAF_file.py
```

```    
python gather_Human_Pfam_Annots.py
```

```    
python gather_Pfam_Clan_Information.py
```

```    
python gather_Combined_RBPs.py
```

```
python plot_Human_SR-SK_range_ProteinFrequencyHeatmap.py
```

```    
python get_SR_prots_and_Plot_RBPproportions_Heatmap.py
```

This series of commands generates Fig 1A, Fig 1B, and Fig S1.
</br></br></br>
For analyses of the features of RS domains (e.g. composition, PTMs, isoforms), run the following commands in-sequence:

```    
python get_PTMs_and_map_to_UniprotIDs.py
```

```    
python calculate_RSdomain_Features.py
```

```    
python plot_RSHKDE_Compositions_Human_RSdomains.py
```

```    
python plot_PTM_and_SRdipep_statistics.py
```

```    
python plot_FullComposition_Human_RSdomains.py
```

This series of commands generates Table S1 and all figure panels appearing in Fig 2, Fig S4, and Fig S5. NOTE: the legend/key is only shown in the Fig 2A file but this legend indicates the categories for all figure panels for these three figures.
</br></br></br>
For comparison of RS domains across human protein isoforms, run the following commands in-sequence:

```
python make_HumanIsoformMap_Uniprot_dict.py
```

```
python get_SR_prots_ALL_ISOFORMS.py
```

```
python analyze_HumanIsoform_RSdomainOverlap.py
```

```
python plot_PieChart_RSdomains_IsoformComparison.py
```

```
python plot_IsoformMatchingProportions.py
```

```
python gather_SelectProtSeqs_IsoformDomainSchematics.py
```

```
python plot_IsoformsDomainSchematics.py
```

This series of commands generates Fig 3, Fig S6, and Table S7.
</br></br></br>
For comparison of human SR/SR-related proteins in this study with those of other studies, run the following commands in-sequence:

```
python plot_VennDiagrams_VariedCompositionThreshold.py
```

```
python plot_DatasetOverlaps_VennDiagrams.py
```

```
python plot_VariedCombinedSRcomposition_LCprots_vs_NewSRprots.py
```

This series of commands generates Table S5, Fig 1C, Fig S3A, Fig S3C, Fig S2B and all venn diagrams appearing in Fig S2A.
</br></br></br>
For in-depth comparison of experimental data for known-, new-, and non-SR-related RBPs:
1. Download "Supplementary information 3(table) from Gerstberger et al. (2014)(PMID:25365966, https://www.nature.com/articles/nrg3813).
2. Open the file in MS Excel, navigate to the spreadsheet titled "RBP table", and save this sheet as the file type "Text (Tab delimited) (\*.txt)" with the file name "SupplementaryInfo3_RNAtarget_Categories.txt".
3. Download the "catrapid_human_basic.zip" file and the "RNAct_supporting_tables.zip" file from the RNAct database (https://rnact.crg.eu/download), then extract the files (NOTE: the "catrapid_human_basic.zip" file will be ~63GB when de-compressed).
4. Run the following commands in-sequence:

```
python plot_RNAtargets_ALL_RBPs.py
```

```
python filter_catRAPID_results_AllProteins_ZSCORE_GREATER_THAN_1_ONLY.py
```

```
python calculate_num_transcripts_bound_per_protein.py
```

This series of commands generates all panels in Fig 4.
