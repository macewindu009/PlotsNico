#include "TFile.h"
#include "TTree.h"


#include <iostream>
#include <fstream>


using namespace std;


int main(int argc, char* argv[] )
{
	cout << "Loading Rootfile.." << endl;
	
	
	
	std::string filedir = argv[1];
	
	cout << "FileDir: " << filedir << endl;

	
	
	std::string longtree = (filedir + "LongZCorrectedRecoil.root").c_str();
	
	TFile *inputFileLong = TFile::Open(longtree.c_str());
	TTree *inputTreeLong = (TTree*)(inputFileLong->Get("t"));
	
	std::string weighttree = (filedir + "flatPtWeight.root").c_str();
	
	TFile *inputFileWeight = TFile::Open(weighttree.c_str());
	TTree *inputTreeWeight = (TTree*)(inputFileWeight->Get("flatPtWeight;2"));
	
	/*
	std::string phitree = (filedir + "PhiCorrectedRecoil.root").c_str();
	
	TFile *inputFilePhi = TFile::Open(phitree.c_str());
	TTree *inputTreePhi = (TTree*)(inputFilePhi->Get("t"));
	*/
	
	
	/*
	cout << "Phi Entries: " << inputTreePhi->GetEntries() << endl;
	cout << "Long Entries: " << inputTreeLong->GetEntries() << endl;
	cout << "Weight Entries: " << inputTreeWeight->GetEntries() << endl;
	if ((inputTreePhi->GetEntries() != inputTreeLong->GetEntries()) && (inputTreeLong->GetEntries() != inputTreeWeight->GetEntries()))
	{
	    cout << "Different amount of entries in tree and Tree with LongZ variables, wrong input data? (quit)";
	    return 1;
	}
	*/
	//TObjArray *entryarrayPhi = new TObjArray(*inputTreePhi->GetListOfBranches());
	
	
	
	
	
	TObjArray *entryarrayLong = new TObjArray(*inputTreeLong->GetListOfBranches());
	
	TObjArray *entryarrayWeight = new TObjArray(*inputTreeWeight->GetListOfBranches());
	
	TObjArray *entryarrayLongOnly = new TObjArray();
	
	//cout << "Number of objects in Phi Tree  = " << entryarrayPhi->GetEntries() << endl; 
	cout << "Number of objects in Long Tree  = " << entryarrayLong->GetEntries() << endl; 
	cout << "Number of objects in Weight Tree  = " << entryarrayWeight->GetEntries() << endl; 
	
	
	for (int i = 0; i < entryarrayLong->GetEntries(); i++)
	    //if (entryarrayPhi->FindObject(entryarrayLong->At(i)->GetName()) == 0)
		entryarrayLongOnly->Add(entryarrayLong->At(i));
	
	
	cout << "Number of objects in Long Only Tree  = " << entryarrayLongOnly->GetEntries() << endl; 
	
	
	TFile *lOFile = new TFile("FullTree.root", "RECREATE");
	TDirectory *cdMAPAnalyzer = lOFile->mkdir("MAPAnalyzer");
	cdMAPAnalyzer->cd();
	
	TTree *OutTree = new TTree();
	
	
	float longTreeEntry[entryarrayLongOnly->GetEntries()];
	for (int i = 0; i < entryarrayLongOnly->GetEntries(); i++)
	{
	    inputTreeLong->SetBranchAddress(entryarrayLongOnly->At(i)->GetName(),&longTreeEntry[i]);
	    OutTree->Branch(entryarrayLongOnly->At(i)->GetName(),&longTreeEntry[i],entryarrayLongOnly->At(i)->GetName()+TString("/F"));
	}
	
	float weightTreeEntry[entryarrayWeight->GetEntries()];
	for (int i = 0; i < entryarrayWeight->GetEntries(); i++)
	{
	    inputTreeWeight->SetBranchAddress(entryarrayWeight->At(i)->GetName(),&weightTreeEntry[i]);
	    OutTree->Branch(entryarrayWeight->At(i)->GetName(),&weightTreeEntry[i],entryarrayWeight->At(i)->GetName()+TString("/F"));
	}
	/*
	float phiTreeEntry[entryarrayPhi->GetEntries()];
	for (int i = 0; i < entryarrayPhi->GetEntries(); i++)
	{
	    inputTreePhi->SetBranchAddress(entryarrayPhi->At(i)->GetName(),&phiTreeEntry[i]);
	    OutTree->Branch(entryarrayPhi->At(i)->GetName(),&phiTreeEntry[i],entryarrayPhi->At(i)->GetName()+TString("/F"));
	}
	*/
	
	
	cout << "Number of objects in OutTree: " << OutTree->GetListOfBranches()->GetEntries() << endl;
	
	int lNEvents = 100;
	for (Long64_t i0=0; i0<lNEvents;i0++)
	{
	    if(i0%1000==0)
		std::cout << "Event " << i0 << "/" << lNEvents << std::endl;
	    //cout << "phi" << endl;
	    //inputTreePhi->GetEntry(i0);
	    //cout << "long" << endl;
	    inputTreeLong->GetEntry(i0);
	    //cout << "weight" << endl;
	    inputTreeWeight->GetEntry(i0);
	    OutTree->Fill();
	}
	
	OutTree->SetName(inputTreeLong->GetName());
	OutTree->Write();
	lOFile->Close();
	
	
	return 0;
}