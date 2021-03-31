# Reproducing SR and SK searches (and subsequent analyses) for the human proteome

### Instructions
1. Download LCD-Composer.py from https://github.com/RossLabCSU/LCD-Composer, as well as all files in the Human_Composition_Searches directory
2. Extract proteome files from compressed file in the same location as LCD-Composer.py
3. Navigate to appropriate folder via command line.
4. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):

```    
python make_Archaea_LCD-Composer_BatchFile.py
```

```
.\RUN_LCD-Composer_Archaea_Batch.bat
```

```
python calculate_Archaea_ProteomeCharacteristics_MeanMedianLCDcontent.py
```

```
python plot_Archaea_Top_LCDcoverage_Proteomes.py
```

This series of commands generates data appearing in Fig 4, Fig 5A, and Table S1.
