inputdir=/storage/jbod/nzaeh/ArtusRun3003NafMVA/merged/
outputdir=/home/nzaeh/MVAPlots/Plots/MT-WJets-Performance/


# MVA MET plots

harry.py -j MVAMET_template.json MVAMET_Transverse_Resolution_pt.json MVAMET_Transverse_Resolution.json MVAMET_resolutionplot.json MVAMET_ptplot.json -d $inputdir -o $outputdir
harry.py -j MVAMET_template.json MVAMET_Transverse_Resolution_NPV.json MVAMET_Transverse_Resolution.json MVAMET_resolutionplot.json MVAMET_npvplot.json -d $inputdir -o $outputdir

harry.py -j MVAMET_template.json MVAMET_Longitudinal_Resolution_pt.json MVAMET_Longitudinal_Resolution.json MVAMET_resolutionplot.json MVAMET_ptplot.json --legend 'lower right' -d $inputdir -o $outputdir
harry.py -j MVAMET_template.json MVAMET_Longitudinal_Resolution_NPV.json MVAMET_Longitudinal_Resolution.json MVAMET_resolutionplot.json MVAMET_npvplot.json -d $inputdir -o $outputdir

