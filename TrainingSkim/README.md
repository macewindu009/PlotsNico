# Using the training data

Once data is skimmed with the SkimManager.py (https://github.com/cms-analysis/HiggsAnalysis-KITHiggsToTauTau/wiki/Performing-a-Skim), it can be processed with tools in this folder. 
Make sure to turn on the usage of MVAMET in the kappa config file.

### Creating plotting files

All createWHATEVERData.sh scripts merge all files of the corresponding decay process, shrink them and produce the corresponding csv file.
They can be run with

```bash
bash makeDYData.sh
```

### Shrinking

shrink.py applies cuts on the skimmed data from SkimManager and saves it in a way that it is able to be processed by makeCSVFromArtus.C


### Make CSV Files

makeCSVFile.C creates a CSV file from data which was produced by the MapAnalyzer (the data which is used for the BDT Training).

makeCSVFromArtus.C creates a CSV file from skimmed data from the SkimManager.py

