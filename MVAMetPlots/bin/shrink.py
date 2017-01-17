#!/usr/bin/env python

from ROOT import *
import sys
import argparse

def Treecopy(folderName, inputfile, outputfile):

    tree = inputfile.Get((folderName + "/ntuple"))

    if tree is None:
        print("No tree found!")
        sys.exit(42)
        
        
    
    branches = ["mvamet","met",
                "recoMetPar","recoMetPerp","recoMetPhi",
                "recoPfMetPar","recoPfMetPerp","recoPfMetPhi",
                "recoilPar","recoilPerp","recoilPhi",
                "pfrecoilPar","pfrecoilPerp","pfrecoilPhi",
                "recoMetOnGenMetProjectionPar","recoMetOnGenMetProjectionPerp","recoMetOnGenMetProjectionPhi",
                "recoPfMetOnGenMetProjectionPar","recoPfMetOnGenMetProjectionPerp","recoPfMetOnGenMetProjectionPhi",    
                "genMetSumEt","genMetPt","genMetPhi",
                "npv","njets","iso_1","iso_2"
                    ]

    for b in tree.GetListOfBranches():
        tree.SetBranchStatus(b.GetName(),0)
        
    for b in branches:
        tree.SetBranchStatus(b,1)

    

    outputfile.cd()
    outputfile.mkdir(folderName)
    outputfile.cd(folderName)

    tree2 = tree.CloneTree()

    tree2.Write()

    return outputfile








if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Shrink Kappa files down to necessary plotting variables of MVAMet.')
    parser.add_argument('-i', '--inputfile', default='/storage/jbod/nzaeh/2016-12-13/Gridoutput/MC.root', help='Path to input root file to be skimmed.')
    parser.add_argument('-o', '--outputfile', default='/storage/jbod/nzaeh/2016-12-13/Gridoutput/MCslimmed.root', help='Path to outputfile/directory.')
    parser.add_argument('-c', '--channels', nargs='*', default=['mm', 'mt'], help='Foldername in which to search for ntuples.')
    args = parser.parse_args()
    parser.add_argument('-m', '--filemode', default='RECREATE', help='create new or add. RECREATE = new file, UPDATE = add to file')
    args = parser.parse_args()

    print('Skimming %s in folder %s in %s mode and saving in %s' %(args.inputfile,args.channels,args.filemode,args.outputfile))

    inputfile = TFile.Open(args.inputfile)
    

    fileout = TFile.Open(args.outputfile,args.filemode)
    
    for channel in args.channels:
        fileout = Treecopy(channel+'_jecUncNom', inputfile, fileout)
    
    fileout.Close()
