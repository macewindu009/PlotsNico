#!/usr/bin/env python

import argparse
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import sys
import ntpath
import os
from array import array

def Treecopy(channel, inputfile, outputfileRoot, expression):

	folderName = channel+'_jecUncNom'
	tree = inputfile.Get((folderName + "/ntuple"))

	if tree is None:
		print("No tree found!")
		sys.exit(42)
	branches = ["mvamet","mvametphi", 
							"met","metphi",
							"mvaMetSumEt", "pfMetSumEt",
							"recoMetPar","recoMetPerp","recoMetPhi",
							"recoPfMetPar","recoPfMetPerp","recoPfMetPhi",
							"recoilPar","recoilPerp","recoilPhi",
							"pfrecoilPar","pfrecoilPerp","pfrecoilPhi",
							"recoMetOnGenMetProjectionPar","recoMetOnGenMetProjectionPerp","recoMetOnGenMetProjectionPhi",
							"recoPfMetOnGenMetProjectionPar","recoPfMetOnGenMetProjectionPerp","recoPfMetOnGenMetProjectionPhi",		
							"genMetSumEt","genMetPt","genMetPhi",
							"npv","npu","njets","iso_1","iso_2","ptvis", "probChiSquare", "probChiSquarePf"
									]

	for b in tree.GetListOfBranches():
		tree.SetBranchStatus(b.GetName(),0)

	for b in branches:
		tree.SetBranchStatus(b,1)


	outputfileRoot.cd()
	outputfileRoot.mkdir(folderName)
	outputfileRoot.cd(folderName)

	tree2 = tree.CloneTree(0)


	for b in tree.GetListOfBranches():
		tree.SetBranchStatus(b.GetName(),1)

	nEntries = tree.GetEntries();

	if expression == '':
		if channel == 'mm':
			expression =	"(extraelec_veto < 0.5)*(iso_1 < 0.15)*(iso_2 < 0.15)*(m_vis > 60.0)*(m_vis < 120.0)"
			#expression =  "(extraelec_veto < 0.5)*(extramuon_veto < 0.5)*(m_vis > 60.0)*(m_vis < 120.0)"
		elif channel == 'mt':
			expression = "(mt_1<40.0)*(againstElectronVLooseMVA6_2 > 0.5)*(againstMuonTight3_2 > 0.5)*(extraelec_veto < 0.5)*(extramuon_veto < 0.5)*(dilepton_veto < 0.5)*(iso_1 < 0.15)*(byTightIsolationMVArun2v1DBoldDMwLT_2 > 0.5)"

	print('Applying weight string: %s'%expression)

	Cut = ROOT.TTreeFormula("name",expression,tree) 
	eventWeight = array('f', [-666.0])
	tree2.Branch('eventWeight',eventWeight, 'eventWeight/F')

	Jet0PtExpr = ROOT.TTreeFormula("name","leadingJetLV.pt()",tree) 
	Jet0_Pt = array('f', [-666.0])
	tree2.Branch('Jet0_Pt',Jet0_Pt, 'Jet0_Pt/F')

	Jet0PhiExpr = ROOT.TTreeFormula("name","leadingJetLV.phi()",tree) 
	Jet0_Phi = array('f', [-666.0])
	tree2.Branch('Jet0_Phi',Jet0_Phi, 'Jet0_Phi/F')

	Jet1PtExpr = ROOT.TTreeFormula("name","trailingJetLV.pt()",tree) 
	Jet1_Pt = array('f', [-666.0])
	tree2.Branch('Jet1_Pt',Jet1_Pt, 'Jet1_Pt/F')

	Jet1PhiExpr = ROOT.TTreeFormula("name","trailingJetLV.phi()",tree) 
	Jet1_Phi = array('f', [-666.0])
	tree2.Branch('Jet1_Phi',Jet1_Phi, 'Jet1_Phi/F')


	for i in range(nEntries):
		tree.GetEntry(i)
		#--- Only write out certain events that pass some cut
		#if passCut(ch.RunNumber, ch.lbn):
		eventWeight[0] = Cut.EvalInstance()
		Jet0_Pt[0] = Jet0PtExpr.EvalInstance()
		Jet0_Phi[0] = Jet0PhiExpr.EvalInstance()
		Jet1_Pt[0] = Jet1PtExpr.EvalInstance()
		Jet1_Phi[0] = Jet1PhiExpr.EvalInstance()
		if eventWeight[0] > 0:
			tree2.Fill()
		if i%100000 == 0:
			print("%i of %i events processed"%(i,nEntries))

	print("New tree in channel %s has %i entries"%(channel,tree2.GetEntries()))

	"""
	if (strcmp( argv[2], "") == 0)
		for (Long64_t i0=0; i0<lNEvents;i0++)
		{
				if(i0%1000000==0)
	std::cout << "Event " << i0 << "/" << lNEvents << std::endl;
				tree->GetEntry(i0);
				tree2->Fill();
		}
	"""
	tree2.Write()

	return outputfileRoot








if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Shrink Kappa files down to necessary plotting variables of MVAMet.')
	parser.add_argument('-i', '--inputfile', default='/storage/jbod/nzaeh/2016-12-13/Gridoutput/merged/MC.root', help='Path to input root file to be skimmed.')
	parser.add_argument('-w', '--weights', default='', help='Weights to apply and save in eventWeights. Default values are standard cuts on mm and mt.')
	parser.add_argument('-p', '--process', default='', help='Use predefined weightstrings per process. Implemented are DY (mm,mt), WJets and Higgs.')
	parser.add_argument('-o', '--outputdir', default='<inputdir>', help='Path to outputdirectory.')
	parser.add_argument('-of', '--outputfile', default='<oldname>_slimmed.root', help='Outputfilename.')
	parser.add_argument('-c', '--channels', nargs='*', default=['mm', 'mt'], help='Foldername in which to search for ntuples.')
	parser.add_argument('-m', '--filemode', default='RECREATE', help='create new or add. RECREATE = new file, UPDATE = add to file')
	args = parser.parse_args()


	if args.outputdir == '<inputdir>':
		outputfiledir = os.path.dirname(args.inputfile) + "/"
	else:
		outputfiledir = args.outputdir
		if not os.path.exists(args.outputdir):
			 os.makedirs(args.outputdir)

	if not args.weights == '':
		expression = args.weights
	else:
		expression = ''



	if args.outputfile == '<oldname>_slimmed.root':
		outputfilename = ntpath.basename(args.inputfile)
		outputfilename = outputfilename[:-5]
		outputfilename = outputfilename + "_slimmed.root"
	else:
		outputfilename = args.outputfile

	outputfile = outputfiledir + outputfilename
	print('Skimming %s in channel %s in %s mode and saving in %s' %(args.inputfile,args.channels,args.filemode,outputfile))


	inputfile = ROOT.TFile.Open(args.inputfile)


	fileout = ROOT.TFile.Open(outputfile,args.filemode)

	for channel in args.channels:
		if not args.process == '':
			if args.process == 'DY':
				if channel == 'mm':
					expression = "(1.0)*(((1.0)))*eventWeight*(((genbosonmass >= 50.0 && (npartons == 0 || npartons >= 5))*3.95423374e-5) + ((genbosonmass >= 50.0 && npartons == 1)*1.27486147e-5) + ((genbosonmass >= 50.0 && npartons == 2)*1.3012785e-5) + ((genbosonmass >= 50.0 && npartons == 3)*1.33802133e-5) + ((genbosonmass >= 50.0 && npartons == 4)*1.09698723e-5)+((genbosonmass < 50.0)*numberGeneratedEventsWeight*crossSectionPerEventWeight))/(numberGeneratedEventsWeight*crossSectionPerEventWeight*sampleStitchingWeight)*(gen_match_1 < 4 || gen_match_2 < 4)*(extraelec_veto < 0.5)*(extramuon_veto < 0.5)*(m_vis > 60.0)*(m_vis < 120.0)*(iso_2 < 0.15)*(iso_1 < 0.15)*((q_1*q_2)<0.0)*zPtReweightWeight"
				elif channel == 'mt':
					expression = "(gen_match_2 == 5)*((((1.0))))*(eventWeight)*((((genbosonmass >= 150.0 && (npartons == 0 || npartons >= 5))*3.95423374e-5) + ((genbosonmass >= 150.0 && npartons == 1)*1.27486147e-5) + ((genbosonmass >= 150.0 && npartons == 2)*1.3012785e-5) + ((genbosonmass >= 150.0 && npartons == 3)*1.33802133e-5) + ((genbosonmass >= 150.0 && npartons == 4)*1.09698723e-5)+((genbosonmass >= 50.0 && genbosonmass < 150.0 && (npartons == 0 || npartons >= 5))*3.95423374e-5) + ((genbosonmass >= 50.0 && genbosonmass < 150.0 && npartons == 1)*1.27486147e-5) + ((genbosonmass >= 50.0 && genbosonmass < 150.0 && npartons == 2)*1.3012785e-5) + ((genbosonmass >= 50.0 && genbosonmass < 150.0 && npartons == 3)*1.33802133e-5) + ((genbosonmass >= 50.0 && genbosonmass < 150.0 && npartons == 4)*1.09698723e-5)+((genbosonmass < 50.0)*numberGeneratedEventsWeight*crossSectionPerEventWeight))/(numberGeneratedEventsWeight*crossSectionPerEventWeight*sampleStitchingWeight)) "
			if args.process == 'WJets':
				if channel == 'mt':
					expression = "(1.0)*(1.0)*(1.0)*(1.0)*(((1.0)*(1.0)))*eventWeight*(((npartons == 0 || npartons >= 5)*7.09390278348407e-4) + ((npartons == 1)*1.90063898596475e-4) + ((npartons == 2)*5.8529964471165e-5) + ((npartons == 3)*1.9206444928444e-5) + ((npartons == 4)*1.923548021385e-5))/(numberGeneratedEventsWeight*crossSectionPerEventWeight*sampleStitchingWeight)*(nbtag == 0)*(extraelec_veto < 0.5)*(extramuon_veto < 0.5)*(againstMuonTight3_2 > 0.5)*(dilepton_veto < 0.5)*(againstElectronVLooseMVA6_2 > 0.5)*(mt_1<40.0)*(byTightIsolationMVArun2v1DBoldDMwLT_2 > 0.5)*(iso_1 < 0.15)*((q_1*q_2)<0.0)"
			if args.process == 'Higgs':
				if channel == 'mt':
					expression = "eventWeight*(nbtag == 0)*(extraelec_veto < 0.5)*(extramuon_veto < 0.5)*(againstMuonTight3_2 > 0.5)*(dilepton_veto < 0.5)*(againstElectronVLooseMVA6_2 > 0.5)*((q_1*q_2)<0.0)"
		else:
			expression = ''
			print('No valid process input. Using no weight cuts')
		fileout = Treecopy(channel, inputfile, fileout, expression)

	fileout.Close()
