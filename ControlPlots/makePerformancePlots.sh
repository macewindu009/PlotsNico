inputdir=/storage/jbod/nzaeh/NAFPF1204NoCorr/merged/
outputdir=/home/nzaeh/MVAPlots/Plots/Uncorr/


# mt wjets performance bash 

folder=mt2016-wjets-performance

harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_pt.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_ptplot.json -d $inputdir -o $outputdir/$folder/
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_NPV.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_npvplot.json -d $inputdir -o $outputdir/$folder/

harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_pt.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_ptplot.json --legend 'lower right' -d $inputdir -o $outputdir/$folder/
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_NPV.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_npvplot.json -d $inputdir -o $outputdir/$folder/

# mt higgs performance bash 

folder=mt2016-higgs-performance

harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_pt.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_ptplot.json -d $inputdir -o $outputdir/$folder/
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_NPV.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_npvplot.json -d $inputdir -o $outputdir/$folder/

harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_pt.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_ptplot.json --legend 'lower right' -d $inputdir -o $outputdir/$folder/
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_NPV.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_npvplot.json -d $inputdir -o $outputdir/$folder/

# mm performance bash 

folder=mm2016-performance

harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_pt.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_ptplot.json -d $inputdir -o $outputdir/$folder/
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Transverse_Resolution_NPV.json $folder/MVAMET_Transverse_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_npvplot.json -d $inputdir -o $outputdir/$folder/

harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_pt.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_ptplot.json --legend 'lower right' -d $inputdir -o $outputdir/$folder/
harry.py -j $folder/MVAMET_template.json $folder/MVAMET_Longitudinal_Resolution_NPV.json $folder/MVAMET_Longitudinal_Resolution.json $folder/MVAMET_resolutionplot.json $folder/MVAMET_npvplot.json -d $inputdir -o $outputdir/$folder/
