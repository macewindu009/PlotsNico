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


bool isUsefulFile(float fileNr, vector<float> filesToInclude)
{
    bool isUseful = false;
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
	TTree *inputTree1 = (TTree*)(inputFile1->Get("t"));

	std::string outputfile = "";

	if (argc > 2)
	    outputfile = argv[2];
	else
	    outputfile = "DataForNNTraining.csv";

	cout << "Outputfile: " <<outputfile << endl;

	fstream dataset;
	dataset.open(outputfile.c_str(), ios::out);


	//define if to filter on filename
	bool applyFilter = false;
	vector<float> fileList = {};
	if (argc > 3)
	{
	    cout << "Files to include: ";
	    for (int i = 3; i < argc; i++)
	    {
		cout << atoi(argv[i]) << ", ";
		fileList.push_back(atoi(argv[i]));
	    }
	    cout << endl;
	    applyFilter = true;
	}
	else
	{
	    cout << "No files specified, not applying any cut on data.";

	    //fileList.push_back(1);
	}

	cout << "Total Number of Files to include: " << fileList.size() << endl;

	cout << "Entries: " << inputTree1->GetEntries() << endl;


	TObjArray *entryarray = new TObjArray(*inputTree1->GetListOfBranches());



	cout << "Number of objects = " << entryarray->GetEntries() << endl;



	cout << "Counting selected events..." << endl;

	int selectedEvents = 0;
	float fileChecker;
	if (applyFilter)
	    inputTree1->SetBranchAddress("fileName",&fileChecker);
	for (int i = 0; i < inputTree1->GetEntries(); i++)
	{
	    inputTree1 -> GetEntry(i);
	    if (!applyFilter)
		selectedEvents++;
	    else
		if (isUsefulFile(fileChecker, fileList))
		    selectedEvents++;

	    if (i % 1000000 == 0)
		cout << i <<  "/" << inputTree1->GetEntries() << " Events processed" << endl;
	}


	int maxEntries = 0;
	bool limitEvents = false;
	if (limitEvents)
	  if (selectedEvents > 500000)
	      maxEntries = 500000;
	  else
	      maxEntries = selectedEvents;
	else
	  maxEntries = selectedEvents;

	cout << "Selected events: " << selectedEvents << endl;

	int variablesize = entryarray->GetEntries();

	cout << "maxEntries: " << maxEntries << endl;
	cout << "variablesize: " << variablesize << endl;
	float *datacontainer = new float[maxEntries*(variablesize)];

	cout << "Loading data into local storage.." << endl;
	int fileID = 0;
	//variablesize - 5 as it does not store the additional 5 calculated target variables
	float varlist[variablesize];
	for (int j = 0; j < entryarray->GetEntries(); j++)
	    {
		inputTree1->SetBranchAddress(entryarray->At(j)->GetName(),&varlist[j]);
		if (applyFilter)
		    if (strcmp(entryarray->At(j)->GetName(),"fileName") == 0)
			fileID = j;

	    }
	int count = 0;
	int treecount = 0;

	while (count < maxEntries && treecount < inputTree1->GetEntries())
	{
	    inputTree1 -> GetEntry(treecount);
	    if (!applyFilter)
	    {
		for (int j = 0; j < variablesize; j++)
		    datacontainer[count*(variablesize)+j] = varlist[j];
		count++;
	    }
	    else
		if (isUsefulFile(varlist[fileID],fileList))
		{
		    for (int j = 0; j < variablesize; j++)
			datacontainer[count*(variablesize)+j] = varlist[j];
		    count++;
		}
	    treecount++;
	    if (treecount % 1000000 == 0)
		cout << treecount << "/" << inputTree1->GetEntries() << " events processed" << endl;
	}
	maxEntries = count;
	cout << "CSV Entries: " << count << endl;
	cout << variablesize << endl;
	float mean[variablesize];
	float variance[variablesize];

	for (int j = 0; j < variablesize; j++)
	{
	    mean[j] = 0;
	    variance[j] = 0;
	    for (int i = 0; i < maxEntries; i++)
		mean[j] += datacontainer[i*(variablesize)+j];
	    mean[j] /= maxEntries;

	    for (int i = 0; i < maxEntries; i++)
		variance[j] += (datacontainer[i*(variablesize)+j] - mean[j])*(datacontainer[i*(variablesize)+j] - mean[j]);
	    variance[j] /= maxEntries;
	}

	int boson_Phi = 0;
	int boson_Pt = 0;
	int recoilslimmedMETs_Phi = 0;
	int recoilslimmedMETs_LongZ = 0;
	int PhiCorrectedRecoil_LongZ = 0;

	for (int i = 0; i < entryarray->GetEntries(); i++)
	{
	    //search for needed variables for target variable
	    if (strcmp(entryarray->At(i)->GetName(),"recoilslimmedMETs_Phi") == 0)
		recoilslimmedMETs_Phi = i;

	    if (strcmp(entryarray->At(i)->GetName(),"recoilslimmedMETs_LongZ") == 0)
		recoilslimmedMETs_LongZ = i;
	    if (strcmp(entryarray->At(i)->GetName(),"Boson_Phi") == 0)
		boson_Phi = i;
	    if (strcmp(entryarray->At(i)->GetName(),"Boson_Pt") == 0)
		boson_Pt = i;
	    if (strcmp(entryarray->At(i)->GetName(),"PhiCorrectedRecoil_LongZ") == 0)
		PhiCorrectedRecoil_LongZ = i;
	    if (strcmp(entryarray->At(i)->GetName(),"fileName") == 0)
	    {
		cout << "FOUND AT LEAST FILENAME" << endl;
		cout << "variance " << variance[i] << endl;
	    }

	    cout << "." << endl;
	    if (variance[i] != 0)
	    {
		if (strcmp(entryarray->At(i)->GetName(),"fileName") == 0)
		    cout << "FOUND FILENAME AND COPYING" << endl;
		dataset <<  entryarray->At(i)->GetName();
		if (i < entryarray->GetEntries()-1)
			dataset << ",";
	    }
	}

	cout << "boson_Phi " << boson_Phi<<endl;
	cout << "boson_Pt " << boson_Pt<<endl;
	cout << "recoilslimmedMETs_Phi " << recoilslimmedMETs_Phi<<endl;
	cout << "recoilslimmedMETs_LongZ " << recoilslimmedMETs_LongZ<<endl;
	cout << "PhiCorrectedRecoil_LongZ " << PhiCorrectedRecoil_LongZ<<endl;

	//add target variables
	dataset << endl;

	cout << "Processing data, generating csv-File.." << endl;
	for (int i = 0; i < maxEntries; i++)
	{
	    for (int j = 0; j < entryarray->GetEntries()-1; j++)
		if (variance[j] != 0)
		    dataset << datacontainer[i*(variablesize)+j] << ",";

	    if (variance[entryarray->GetEntries()-1] != 0)
	        dataset << datacontainer[i*variablesize + entryarray->GetEntries()-1];

	    if (i != (inputTree1->GetEntries()-1))
	    {
		dataset << endl;
	    }

	    if (i % 100000 == 0)
		cout << i <<  "/" << maxEntries << " Events processed" << endl;
	}

	cout << "Succesfully created csv File!" << endl;
	dataset.close();
	return 0;
}
