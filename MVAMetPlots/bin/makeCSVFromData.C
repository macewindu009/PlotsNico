#include "TFile.h"
#include "TTree.h"


#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>


using namespace std;

const float  PI_F=3.14159265358979f;

/* Files belong to the following variable refered to in fileName
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 --> 1
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 --> 8
DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 --> 9
DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 --> 6
DYJetsToLL_M-150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 --> 7
DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 --> 2
DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 --> 3
DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 --> 4
DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 --> 5
*/
/*

bool isUsefulFile(float fileNr)
{
    bool isUseful = false;
    vector<float> filesToInclude = {8};
    for (unsigned int i = 0; i < filesToInclude.size(); i++)
        if (fileNr == filesToInclude[i])
	    isUseful = true;
	
    return isUseful;
  
}
*/

int main(int argc, char* argv[] )
{
	cout << "Loading Rootfile.." << endl;
	
        cout << "argc " << argc << endl;
        for (int i = 3; i < argc;i++)
        {
            cout << "Creating CSV from file " << argv[1] << " in channel " << argv[i] << endl;
            string channelpostfix = "_jecUncNom/ntuple";
            string channelname = argv[i] + channelpostfix;
            cout << channelname << endl;
            std::string filename = argv[1];
            
            TFile *inputFile = TFile::Open(filename.c_str());
            TTree *inputTree = (TTree*)(inputFile->Get(channelname.c_str()));
            cout << "Entries: " << inputTree->GetEntries() << endl;
            
            TObjArray *entryarray = new TObjArray(*inputTree->GetListOfBranches());
            cout << "Number of objects = " << entryarray->GetEntries() << endl; 
            
            int beginIdx = filename.rfind('/');
            string filenameout = filename.substr(beginIdx + 1);
            filenameout.erase (filenameout.end()-5, filenameout.end());
            string channel = argv[i];
            string outputdir = argv[2];
            filenameout = outputdir + "/" + channel + "_" + filenameout + ".csv";
            cout << "CSV-filename: " << filenameout << endl;
            fstream dataset;
            
            
            
            //string fileoutname = filename + ""
            dataset.open(filenameout, ios::out);
            
            int maxEntries = 0;
            bool limitEvents = false;
            if (limitEvents)
            if (inputTree->GetEntries() > 100000)
                maxEntries = 100000;
            else 
                maxEntries = inputTree->GetEntries();
            else
                maxEntries = inputTree->GetEntries();
            //cout << maxEntries << endl;
            
            vector<string> branches = {"mvamet","met",
                "recoMetPar","recoMetPerp","recoMetPhi",
                "recoPfMetPar","recoPfMetPerp","recoPfMetPhi",
                "recoilPar","recoilPerp","recoilPhi",
                "pfrecoilPar","pfrecoilPerp","pfrecoilPhi",
                "recoMetOnGenMetProjectionPar","recoMetOnGenMetProjectionPerp","recoMetOnGenMetProjectionPhi",
                "recoPfMetOnGenMetProjectionPar","recoPfMetOnGenMetProjectionPerp","recoPfMetOnGenMetProjectionPhi",    
                "genMetSumEt","genMetPt","genMetPhi",
                "npv","njets","iso_1","iso_2", "ptvis"
            };
                    
            /*
            for (unsigned int j = 0; j < branches.size(); j++)
                cout << branches[j] << endl;
            */
            int variablesize = branches.size();
	
            vector<float> datacontainer(maxEntries*(variablesize));
            
            cout << "Loading data into local storage.." << endl;

            float varlist[(variablesize)];
            int varlistInt[variablesize];
            int typelist[variablesize];
            for (int j = 0; j < variablesize; j++)
                {
                    typelist[j] = inputTree->SetBranchAddress(branches[j].c_str(),&varlist[j]);
                    //cout << inputTree->SetBranchAddress(branches[j].c_str(),&varlist[j]) << endl;
                    if (typelist[j] != 0)
                    {
                        inputTree->SetBranchAddress(branches[j].c_str(),&varlistInt[j]);
                        cout << "Reset Branch address for integer values, should be fixed now." << endl;
                    }
                }
            
            int count = 0;
            //for (int i = 0; i < maxEntries; i++)
            while (count < maxEntries)
            {
                inputTree -> GetEntry(count);
                for (int j = 0; j < variablesize; j++)
                {
                    if (typelist[j] == 0)
                        datacontainer[count*(variablesize)+j] = varlist[j];
                    else
                        datacontainer[count*(variablesize)+j] = float(varlistInt[j]);
                    
                }
                count++;
                if (count % 100000 == 0)
                    cout << count << " of " << maxEntries << " Events processed" << endl;
            }	
            maxEntries = count;
            cout << "CSV Entries: " << count << endl;
            
            for (int j = 0; j < variablesize; j++)
            {
                dataset <<  branches[j] << ",";
            }
            
            dataset << endl;
	
            cout << "Processing data, generating csv-File.." << endl;
            for (int j = 0; j < maxEntries; j++)
            {
                for (int k = 0; k < variablesize-1; k++)
                    dataset << datacontainer[j*(variablesize)+k] << ",";
                dataset << datacontainer[j*(variablesize)+variablesize-1];
                
                if (j != (maxEntries-1))
                {
                    dataset << endl;
                }
                
                if (j % 100000 == 0)
                    cout << j << " of " << maxEntries << " Events processed" << endl;
            }
            
            
            dataset.close();
        }    
            
	return 0;
}
