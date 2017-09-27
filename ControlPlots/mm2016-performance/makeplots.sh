inputdir=/storage/jbod/nzaeh/NAFOwnCorr2605/merged/
outputdir=/home/nzaeh/MVAPlots/Plots/NAFOwnCorr2605/


# MVA MET plots

harry.py -j MVAMET_template.json MVAMET_Transverse_Resolution_pt.json MVAMET_Transverse_Resolution.json MVAMET_resolutionplot.json MVAMET_ptplot.json -d $inputdir -o $outputdir --formats 'pdf' & 
harry.py -j MVAMET_template.json MVAMET_Transverse_Resolution_NPV.json MVAMET_Transverse_Resolution.json MVAMET_resolutionplot.json MVAMET_npvplot.json -d $inputdir -o $outputdir --formats 'pdf' &
harry.py -j MVAMET_template.json MVAMET_Transverse_Resolution_NJets.json MVAMET_Transverse_Resolution.json MVAMET_resolutionplot.json MVAMET_NJetsplot.json -d $inputdir -o $outputdir --formats 'pdf' --legend 'upper left' --x-ticks 0 1 2 3 &

harry.py -j MVAMET_template.json MVAMET_Longitudinal_Resolution_pt.json MVAMET_Longitudinal_Resolution.json MVAMET_resolutionplot.json MVAMET_ptplot.json --legend 'lower right' -d $inputdir -o $outputdir --formats 'pdf' &
harry.py -j MVAMET_template.json MVAMET_Longitudinal_Resolution_NPV.json MVAMET_Longitudinal_Resolution.json MVAMET_resolutionplot.json MVAMET_npvplot.json -d $inputdir -o $outputdir --formats 'pdf' &
harry.py -j MVAMET_template.json MVAMET_Longitudinal_Resolution_NJets.json MVAMET_Longitudinal_Resolution.json MVAMET_resolutionplot.json MVAMET_NJetsplot.json -d $inputdir -o $outputdir --formats 'pdf' --legend 'upper left' --x-ticks 0 1 2 3

