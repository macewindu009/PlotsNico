mergedDir=/storage/jbod/nzaeh/NAFMVA1204NoCorr/merged
saveDir=$mergedDir/PlotData/
cd $mergedDir
if [ ! -d "$saveDir" ]; then
# Control will enter here if $saveDir doesn't exist.
	mkdir $saveDir
fi
hadd -f $saveDir/WJets.root W1JetsToLNu_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8/*.root W2JetsToLNu_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8/*.root W2JetsToLNu_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8_ext1/*.root W3JetsToLNu_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8/*.root W3JetsToLNu_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8_ext1/*.root W4JetsToLNu_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8/*.root W4JetsToLNu_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8_ext1/*.root W4JetsToLNu_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8_ext2/*.root WJetsToLNu_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8/*.root WJetsToLNu_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8_ext2/*.root
cd -
python shrink.py -i $saveDir/WJets.root -o $saveDir -c mt -p WJets
makeCSVFromArtus $saveDir/WJets_slimmed.root $saveDir mt
