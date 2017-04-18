# MVA MET plots
harry.py -j MVAMET_Response.json MVAMET_Response_pt.json MVAMET_template.json MVAMET_ptplot.json
harry.py -j MVAMET_Response.json MVAMET_Response_NPV.json MVAMET_template.json MVAMET_npvplot.json

harry.py -j MVAMET_template.json MVAMET_Transverse_Resolution_pt.json MVAMET_Transverse_Resolution.json MVAMET_resolutionplot.json MVAMET_ptplot.json --nicks-blacklist function
harry.py -j MVAMET_resolutionplot_fit.json MVAMET_template.json MVAMET_Transverse_Resolution_NPV.json MVAMET_Transverse_Resolution.json  MVAMET_npvplot.json 

harry.py -j MVAMET_template.json MVAMET_Longitudinal_Resolution_pt.json MVAMET_Longitudinal_Resolution.json MVAMET_resolutionplot.json MVAMET_ptplot.json --legend 'lower right' --nicks-blacklist function
harry.py -j MVAMET_resolutionplot_fit.json MVAMET_template.json MVAMET_Longitudinal_Resolution_NPV.json MVAMET_Longitudinal_Resolution.json  MVAMET_npvplot.json

# input plots
harry.py -j INPUT_Response.json INPUT_Response_pt.json INPUT_template.json MVAMET_ptplot.json
harry.py -j INPUT_Response.json INPUT_Response_NPV.json INPUT_template.json MVAMET_npvplot.json

harry.py -j INPUT_template.json INPUT_Transverse_Resolution_pt.json INPUT_Transverse_Resolution.json INPUT_resolutionplot.json MVAMET_ptplot.json
harry.py -j INPUT_template.json INPUT_Transverse_Resolution_NPV.json INPUT_Transverse_Resolution.json INPUT_resolutionplot.json MVAMET_npvplot.json

harry.py -j INPUT_template.json INPUT_Longitudinal_Resolution_pt.json INPUT_Longitudinal_Resolution.json INPUT_resolutionplot.json MVAMET_ptplot.json --legend 'lower right'
harry.py -j INPUT_template.json INPUT_Longitudinal_Resolution_NPV.json INPUT_Longitudinal_Resolution.json INPUT_resolutionplot.json MVAMET_npvplot.json
