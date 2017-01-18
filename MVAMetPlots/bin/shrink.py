#!/usr/bin/env python

import argparse
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import sys
import ntpath
import os

def Treecopy(channel, inputfile, outputfileRoot):

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
                "npv","npu","njets","iso_1","iso_2","ptvis"
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
    if channel == 'mm':
        expression =  "(extraelec_veto < 0.5)*(iso_1 < 0.15)*(iso_2 < 0.15)*(m_vis > 60.0)*(m_vis < 120.0)"
        #expression =  "(extraelec_veto < 0.5)*(extramuon_veto < 0.5)*(m_vis > 60.0)*(m_vis < 120.0)"
    elif channel == 'mt':
        expression = "(mt_1<40.0)*(againstElectronVLooseMVA6_2 > 0.5)*(againstMuonTight3_2 > 0.5)*(extraelec_veto < 0.5)*(extramuon_veto < 0.5)*(dilepton_veto < 0.5)*(iso_1 < 0.15)*(byTightIsolationMVArun2v1DBoldDMwLT_2 > 0.5)"
    else:
        expression = ""

            
    
    Cut = ROOT.TTreeFormula("name",expression,tree) 
    for i in range(nEntries):
        tree.GetEntry(i)
        #--- Only write out certain events that pass some cut
        #if passCut(ch.RunNumber, ch.lbn):
	result = Cut.EvalInstance()
	if result > 0:
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
    parser.add_argument('-o', '--outputdir', default='<inputdir>', help='Path to outputdirectory.')
    parser.add_argument('-of', '--outputfile', default='<oldname>_slimmed.root', help='Outputfilename.')
    parser.add_argument('-c', '--channels', nargs='*', default=['mm', 'mt'], help='Foldername in which to search for ntuples.')
    parser.add_argument('-m', '--filemode', default='RECREATE', help='create new or add. RECREATE = new file, UPDATE = add to file')
    args = parser.parse_args()

    
    if args.outputdir == '<inputdir>':
        outputfiledir = os.path.dirname(args.inputfile) + "/"
    else:
        outputfiledir = args.outputdir
        
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
        fileout = Treecopy(channel, inputfile, fileout)
    
    fileout.Close()
