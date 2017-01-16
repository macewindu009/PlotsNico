#include "TFile.h"
#include "TTree.h"


#include <iostream>
#include <fstream>
#include <math.h>



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

bool isUsefulFile(float fileNr)
{
    bool isUseful = false;
    vector<float> filesToInclude = {8};
    for (unsigned int i = 0; i < filesToInclude.size(); i++)
        if (fileNr == filesToInclude[i])
	    isUseful = true;
	
    return isUseful;
  
}


int main(int argc, char* argv[] )
{
	cout << "Loading Rootfile.." << endl;
	
	std::string filename1 = argv[1];
	TFile *inputFile1 = TFile::Open(filename1.c_str());
	TTree *inputTree1 = (TTree*)(inputFile1->Get("MAPAnalyzer/t"));
	
	std::string filenameFriend = argv[2];
	
	std::string filenameFriendPhi = filenameFriend + "PhiCorrectedRecoil.root";
	TFile *inputFileFriendPhi = TFile::Open(filenameFriendPhi.c_str());
	TTree *inputTreeFriendPhi = (TTree*)(inputFileFriendPhi->Get("PhiCorrectedRecoil"));

	std::string filenameFriendRecoil = filenameFriend + "LongZCorrectedRecoil.root";
	TFile *inputFileFriendLongZ = TFile::Open(filenameFriendRecoil.c_str());
	TTree *inputTreeFriendLongZ = (TTree*)(inputFileFriendLongZ->Get("LongZCorrectedRecoil"));
	

	
	cout << "Entries: " << inputTree1->GetEntries() << endl;
	cout << "Entries friend Phi: " << inputTreeFriendPhi->GetEntries() << endl;
	cout << "Entries friend LongZ: " << inputTreeFriendLongZ->GetEntries() << endl;
	
	if (inputTree1->GetEntries() != inputTreeFriendPhi->GetEntries() || inputTreeFriendPhi->GetEntries() != inputTreeFriendLongZ->GetEntries())
	{
	    cout << "Different amount of entries in tree and its friends, wrong input data? (quit)";
	    return 1;
	}
	
	TObjArray *entryarray = new TObjArray(*inputTree1->GetListOfBranches());

	TObjArray *entryarrayPhi = new TObjArray(*inputTreeFriendPhi->GetListOfBranches());

	TObjArray *entryarrayLongZ = new TObjArray(*inputTreeFriendLongZ->GetListOfBranches());

	cout << "Number of objects = " << entryarray->GetEntries() << endl; 
	cout << "Number of objects friend Phi = " << entryarrayPhi->GetEntries() << endl; 
	cout << "Number of objects friend LongZ= " << entryarrayLongZ->GetEntries() << endl; 
	//cout << "Anzahl Objekte = " << inputTree1->GetListOfBranches()->GetEntries() << endl; 

	fstream dataset;
	dataset.open("dataMVAMet.csv", ios::out);
	
	

	int maxEntries = 0;
	bool limitEvents = false;
	if (limitEvents)
	  if (inputTree1->GetEntries() > 1000000)
	      maxEntries = 1000000;
	  else 
	      maxEntries = inputTree1->GetEntries();
	else
	  maxEntries = inputTree1->GetEntries();

	int variablesize = entryarray->GetEntries()+3+entryarrayPhi->GetEntries()+entryarrayLongZ->GetEntries();
	
	vector<float> datacontainer(maxEntries*(variablesize));
	
	cout << "Loading data into local storage.." << endl;
	int selectentry = 0;
	//variablesize - 3 as it does not store the additional 3 calculated target variables
	float varlist[(variablesize-3)];
	for (int j = 0; j < entryarray->GetEntries(); j++)
	    {
		inputTree1->SetBranchAddress(entryarray->At(j)->GetName(),&varlist[j]);
		if (strcmp(entryarray->At(j)->GetName(),"select") == 0)
		    selectentry = j;
	    }
	    
	for (int j = 0; j < entryarrayPhi->GetEntries(); j++)
	    {
		inputTreeFriendPhi->SetBranchAddress(entryarrayPhi->At(j)->GetName(),&varlist[j+entryarray->GetEntries()]);
	    }
	    
	for (int j = 0; j < entryarrayLongZ->GetEntries(); j++)
	    {
		inputTreeFriendLongZ->SetBranchAddress(entryarrayLongZ->At(j)->GetName(),&varlist[j+entryarray->GetEntries()+entryarrayPhi->GetEntries()]);
	    }
	    
	
	int count = 0;
	int treecount = 0;
	//for (int i = 0; i < maxEntries; i++)
	while (count < maxEntries && treecount < inputTree1->GetEntries())
	{
	    inputTree1 -> GetEntry(treecount);
	    inputTreeFriendLongZ -> GetEntry(treecount);
	    inputTreeFriendPhi -> GetEntry(treecount);
	    
	    if (varlist[selectentry] == 1)
	    {
		for (int j = 0; j < variablesize-3; j++)
		{
		    
		    {
			datacontainer[count*(variablesize)+j] = varlist[j];	
		    }
		}
		count++;
	    }
	    treecount++;
	}	
	maxEntries = count;
	cout << "CSV Entries: " << count << endl;
	
	float mean[variablesize-3];
	float variance[variablesize-3];
	
	for (int j = 0; j < variablesize-3; j++)
	{
	    mean[j] = 0;
	    variance[j] = 0;
	    for (int i = 0; i < maxEntries; i++)
		mean[j] += datacontainer[i*(variablesize)+j];
	    mean[j] /= maxEntries;
	    
	    for (int i = 0; i < maxEntries; i++)
		variance[j] = (datacontainer[i*(variablesize)+j] - mean[j])*(datacontainer[i*(variablesize)+j] - mean[j]);
	    variance[j] /= maxEntries;
	}

	int boson_Phi = 0;
	int boson_Pt = 0;
	int recoilslimmedMETs_Phi = 0;
	int recoilslimmedMETs_Pt = 0;
	int PhiCorrectedRecoil_LongZ = 0;
	
	for (int i = 0; i < entryarray->GetEntries(); i++)
	{
	    //search for needed variables for target variable
	    if (strcmp(entryarray->At(i)->GetName(),"recoilslimmedMETs_Phi") == 0)
		recoilslimmedMETs_Phi = i;
	    if (strcmp(entryarray->At(i)->GetName(),"recoilslimmedMETs_Pt") == 0)
		recoilslimmedMETs_Phi = i;
	    if (strcmp(entryarray->At(i)->GetName(),"Boson_Phi") == 0)
		boson_Phi = i;
	    if (strcmp(entryarray->At(i)->GetName(),"Boson_Pt") == 0)
		boson_Pt = i;
	    if (variance[i] != 0)
		dataset <<  entryarray->At(i)->GetName() << ",";
	}
	
	
	
	
	for (int i = 0; i < entryarrayPhi->GetEntries(); i++)
	{
	    if (strcmp(entryarrayPhi->At(i)->GetName(),"PhiCorrectedRecoil_LongZ") == 0)
		PhiCorrectedRecoil_LongZ = i;
	    if (variance[i + entryarray->GetEntries()] != 0)
		dataset <<  entryarrayPhi->At(i)->GetName() << ",";
	}
	
	cout << "boson_Phi " << boson_Phi<<endl;
	cout << "boson_Pt " << boson_Pt<<endl;
	cout << "recoilslimmedMETs_Phi " << recoilslimmedMETs_Phi<<endl;
	cout << "recoilslimmedMETs_Pt " << recoilslimmedMETs_Pt<<endl;
	cout << "PhiCorrectedRecoil_LongZ " << PhiCorrectedRecoil_LongZ<<endl;
	
	for (int i = 0; i < entryarrayLongZ->GetEntries(); i++)
	{
	    if (variance[i + entryarray->GetEntries()+entryarrayPhi->GetEntries()] != 0)
		dataset <<  entryarrayLongZ->At(i)->GetName() << ",";
	}
	
	//add target variables
	dataset << "targetPhiFromSlimmed,targetRecoilFromBDT,targetRecoilFromSlimmed";
	dataset << endl;
	
	cout << "Processing data, generating csv-File.." << endl;
	for (int i = 0; i < maxEntries; i++)
	{
	    for (int j = 0; j < entryarray->GetEntries(); j++)
		if (variance[j] != 0)
		    dataset << datacontainer[i*(variablesize)+j] << ",";
	    
	    for (int j = 0; j < entryarrayPhi->GetEntries(); j++)
		if (variance[j+entryarray->GetEntries()] != 0)
		    dataset << datacontainer[i*(variablesize)+j+entryarray->GetEntries()] << ",";
		
	    for (int j = 0; j < entryarrayLongZ->GetEntries(); j++)
		if (variance[j+entryarray->GetEntries()+entryarrayPhi->GetEntries()] != 0)
		    dataset << datacontainer[i*(variablesize)+j+entryarray->GetEntries()+entryarrayPhi->GetEntries()] << ",";
	    //extra first target variable
	    //"(Boson_Phi-recoilslimmedMETs_Phi+TMath::Pi()) - 2.*TMath::Pi()*((Boson_Phi-recoilslimmedMETs_Phi) > 0)",
	    if (datacontainer[i*(variablesize)+boson_Phi]-datacontainer[i*(variablesize)+recoilslimmedMETs_Phi] > 0)
		dataset << datacontainer[i*(variablesize)+boson_Phi]-datacontainer[i*(variablesize)+recoilslimmedMETs_Phi]-PI_F;
	    else
		dataset << datacontainer[i*(variablesize)+boson_Phi]-datacontainer[i*(variablesize)+recoilslimmedMETs_Phi]+PI_F;
	    
	    
	    
	    
	    dataset << ",";
	    //extra second target variable depending on the BDT prediction
	    //"-Boson_Pt/PhiCorrectedRecoil_LongZ",
	    dataset << (-datacontainer[i*(variablesize)+boson_Pt]/datacontainer[i*(variablesize)+entryarray->GetEntries()+PhiCorrectedRecoil_LongZ]);
	    
	    //extra third target Pt variable independent of the BDT prediction
	    //"-Boson_Pt/PhiCorrectedRecoil_LongZ",
	    dataset << (-datacontainer[i*(variablesize)+boson_Pt]/datacontainer[i*(variablesize)+recoilslimmedMETs_Pt]);
	    
	    
	    if (i != (inputTree1->GetEntries()-1))
	    {
		dataset << endl;
	    }
	    
	    if (i % 10000 == 0)
		cout << i << " Events processed" << endl;
	}
	
	
	dataset.close();
	  
	return 0;
}
