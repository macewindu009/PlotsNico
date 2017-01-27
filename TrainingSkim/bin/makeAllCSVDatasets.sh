#!/bin/bash
for trainingname in $( ls ../ConfignamesJet0/ | grep With); do	
	makeNNDatasetFromDataAndWeight /storage/jbod/nzaeh/2016-10-04/Gridoutput/$trainingname/data1.root /storage/jbod/nzaeh/2016-10-04/Gridoutput/$trainingname/flatPtWeight.root /storage/jbod/nzaeh/2016-10-04/Gridoutput/$trainingname/dataMVAMet$trainingname.csv 8
done
