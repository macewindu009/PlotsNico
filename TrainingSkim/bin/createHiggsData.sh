mergedDir=/storage/jbod/nzaeh/NAFMVA1204NoCorr/merged
saveDir=$mergedDir/PlotData/
if [ ! -d "$saveDir" ]; then
# Control will enter here if $saveDir doesn't exist.
	mkdir $saveDir
fi
python shrink.py -i $mergedDir/VBFHToTauTauM125_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_powheg-pythia8/*.root -o $saveDir -of qqHiggs_slimmed.root -c mt -p Higgs
makeCSVFromArtus $saveDir/qqHiggs_slimmed.root $saveDir mt
python shrink.py -i $mergedDir/GluGluHToTauTauM125_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_powheg-pythia8/*.root -o $saveDir -of ggHiggs_slimmed.root -c mt -p Higgs
makeCSVFromArtus $saveDir/ggHiggs_slimmed.root $saveDir mt
