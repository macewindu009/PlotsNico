inputdir=/storage/jbod/nzaeh/NAFMVAOwnCorr2604/merged/
outputdir=/home/nzaeh/MVAPlots/Plots/OwnCorr/
loglevel=""

# mt wjets performance bash 

folder=mt2016-wjets-performance

harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_pt.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_ptplot.json -d $inputdir -o $outputdir/$folder/ $loglevel
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_NPV.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_npvplot.json -d $inputdir -o $outputdir/$folder/ $loglevel
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_NJets.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_NJetsplot.json -d $inputdir -o $outputdir/$folder/ $loglevel

harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_pt.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_ptplot.json --legend 'lower right' -d $inputdir -o $outputdir/$folder/ $loglevel
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_NPV.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_npvplot.json -d $inputdir -o $outputdir/$folder/ $loglevel
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_NJets.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_NJetsplot.json -d $inputdir -o $outputdir/$folder/ $loglevel
# mt higgs performance bash 

folder=mt2016-higgs-performance

harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_pt.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_ptplot.json -d $inputdir -o $outputdir/$folder/ $loglevel
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_NPV.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_npvplot.json -d $inputdir -o $outputdir/$folder/ $loglevel 
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_NJets.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_NJetsplot.json -d $inputdir -o $outputdir/$folder/ $loglevel

harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_pt.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_ptplot.json --legend 'lower right' -d $inputdir -o $outputdir/$folder/ $loglevel
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_NPV.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_npvplot.json -d $inputdir -o $outputdir/$folder/ $loglevel
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_NJets.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_NJetsplot.json -d $inputdir -o $outputdir/$folder/ $loglevel

# mt dy performance bash 

folder=mt2016-DY-performance

harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_pt.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_ptplot.json -d $inputdir -o $outputdir/$folder/ $loglevel
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_NPV.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_npvplot.json -d $inputdir -o $outputdir/$folder/ $loglevel 
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_NJets.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_NJetsplot.json -d $inputdir -o $outputdir/$folder/ $loglevel

harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_pt.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_ptplot.json --legend 'lower right' -d $inputdir -o $outputdir/$folder/ $loglevel
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_NPV.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_npvplot.json -d $inputdir -o $outputdir/$folder/ $loglevel
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_NJets.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_NJetsplot.json -d $inputdir -o $outputdir/$folder/ $loglevel


# mm performance bash 

folder=mm2016-performance

harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_pt.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_ptplot.json -d $inputdir -o $outputdir/$folder/ $loglevel
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_NPV.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_npvplot.json -d $inputdir -o $outputdir/$folder/ $loglevel
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_NJets.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_NJetsplot.json -d $inputdir -o $outputdir/$folder/ $loglevel

harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_pt.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_ptplot.json --legend 'lower right' -d $inputdir -o $outputdir/$folder/ $loglevel
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_NPV.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_npvplot.json -d $inputdir -o $outputdir/$folder/ $loglevel
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_NJets.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_NJetsplot.json -d $inputdir -o $outputdir/$folder/ $loglevel
