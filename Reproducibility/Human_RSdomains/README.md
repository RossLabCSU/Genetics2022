# Reproducing SR and SK searches (and subsequent analyses) for the human proteome

### Instructions
1. Download LCD-Composer.py from https://github.com/RossLabCSU/LCD-Composer, as well as all files in the Human_Composition_Searches directory.
2. Extract all compressed files in the same location as LCD-Composer.py and other Python scripts.
3. Navigate to appropriate folder via command line.
4. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):

```    
python make_Human_SRsearch_BatchFile.py
```

```
.\RUN_LCD-Composer_Human_SRcompRANGE_Batch.bat
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
python get_SR_prots_and_Plot_RBPproportions_Heatmap.py
```

This series of commands generates data appearing in Fig 1A,B and Fig S1.

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
