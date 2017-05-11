mergedDir=/storage/jbod/nzaeh/NAFAlexei0805/merged
saveDir=$mergedDir/PlotData/
cd $mergedDir
if [ ! -d "$saveDir" ]; then
# Control will enter here if $saveDir doesn't exist.
	mkdir $saveDir
fi
hadd -f $saveDir/DY.root DY1JetsToLLM50_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8/*.root DY3JetsToLLM50_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8/*.root DY4JetsToLLM50_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8/*.root DY2JetsToLLM50_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8/*.root DYJetsToLLM10to50_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8/*.root DYJetsToLLM50_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8_ext1/*.root DYJetsToLLM50_RunIISummer16MiniAODv2_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8_ext2/*.root
cd -
#python shrink.py -i $saveDir/DY.root -o $saveDir -c mm mt -p DY
makeCSVFromArtus $saveDir/DY_slimmed.root $saveDir mm mt
