inputdir=/storage/jbod/nzaeh/ArtusRun3003NafMVA/merged/
outputdir=/home/nzaeh/MVAPlots/Plots/MM-PerformanceRogerMean


# MVA MET plots

harry.py -j MVAMET_template.json MVAMET_Transverse_Resolution_pt.json MVAMET_Transverse_Resolution.json MVAMET_resolutionplot.json MVAMET_ptplotRoger0Jet.json -d $inputdir -o $outputdir --filename "RogerTrans0Jet"
harry.py -j MVAMET_template.json MVAMET_Transverse_Resolution_pt.json MVAMET_Transverse_Resolution.json MVAMET_resolutionplot.json MVAMET_ptplotRoger1Jet.json -d $inputdir -o $outputdir --filename "RogerTrans1Jet"

harry.py -j MVAMET_template.json MVAMET_Transverse_Resolution_pt.json MVAMET_Transverse_Resolution.json MVAMET_resolutionplot.json MVAMET_ptplotRoger2Jet.json -d $inputdir -o $outputdir --filename "RogerTrans2Jet"


harry.py -j MVAMET_template.json MVAMET_Longitudinal_Resolution_pt.json MVAMET_Longitudinal_Resolution.json MVAMET_resolutionplot.json MVAMET_ptplotRoger0Jet.json --legend 'lower right' -d $inputdir -o $outputdir --filename "RogerLong0Jet"

harry.py -j MVAMET_template.json MVAMET_Longitudinal_Resolution_pt.json MVAMET_Longitudinal_Resolution.json MVAMET_resolutionplot.json MVAMET_ptplotRoger1Jet.json --legend 'lower right' -d $inputdir -o $outputdir --filename "RogerLong1Jet"

harry.py -j MVAMET_template.json MVAMET_Longitudinal_Resolution_pt.json MVAMET_Longitudinal_Resolution.json MVAMET_resolutionplot.json MVAMET_ptplotRoger2Jet.json --legend 'lower right' -d $inputdir -o $outputdir --filename "RogerLong2Jet"


