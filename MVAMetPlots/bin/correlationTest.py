#!/usr/bin/env python


from __future__ import print_function

import sys
import os
import os.path
import time
#import ROOT

import time
import csv
import numpy as np
import math
import argparse

import json

import matplotlib.mlab as mlab
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm
from matplotlib.mlab import PCA
from scipy.optimize import curve_fit
from scipy.stats import chisquare


"""
Files are referred to in variable fileName as
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 --> 1
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 --> 8
DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 --> 9
DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 --> 6
DYJetsToLL_M-150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 --> 7
DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 --> 2
DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 --> 3
DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 --> 4
DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 --> 5
"""

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def doublegauss(x, *p):
    A1, A2, mu1, mu2, sigma1, sigma2 = p
    return A1*np.exp(-(x-mu1)**2/(2.*sigma1**2)) + A2*np.exp(-(x-mu2)**2/(2.*sigma2**2))


def load_datasetcsv(config):
    #create Treevariable
    
    start = time.time()
    

    ArtusDict = {
            "mvamet" : "dmvamet_Pt",
            "mvametphi" : "dmvamet_Phi",
            "met" : "dpfmet_Pt",
            "metphi" : "dpfmet_Phi",
            "mvaMetSumEt" : "LongZCorrectedRecoil_sumEt",
            "pfMetSumEt" : "recoilslimmedMETs_sumEt",
            "recoMetPar" : "recoMetPar",
            "recoMetPerp" : "recoMetPerp",
            "recoMetPhi" : "recoMetPhi",
            "recoPfMetPar" : "recoPfMetPar",
            "recoPfMetPerp" : "recoPfMetPerp",
            "recoPfMetPhi" : "recoPfMetPhi",
            "recoilPar" : "LongZCorrectedRecoil_LongZ",
            "recoilPerp" : "LongZCorrectedRecoil_PerpZ",
            "recoilPhi" : "LongZCorrectedRecoil_Phi",
            "pfrecoilPar" : "recoilslimmedMETs_LongZ",
            "pfrecoilPerp" : "recoilslimmedMETs_PerpZ",
            "pfrecoilPhi" : "recoilslimmedMETs_Phi",
            "recoMetOnGenMetProjectionPar" : "recoMetOnGenMetProjectionPar",
            "recoMetOnGenMetProjectionPerp" : "recoMetOnGenMetProjectionPerp",
            "recoMetOnGenMetProjectionPhi" : "recoMetOnGenMetProjectionPhi",
            "recoPfMetOnGenMetProjectionPar" : "recoPfMetOnGenMetProjectionPar",
            "recoPfMetOnGenMetProjectionPerp" : "recoPfMetOnGenMetProjectionPerp",
            "recoPfMetOnGenMetProjectionPhi" : "recoPfMetOnGenMetProjectionPhi",     
            "genMetSumEt" : "genMetSumEt",
            "genMetPt" : "genMet_Pt",
            "genMetPhi" : "genMetPhi",
            "npv" : "NVertex",
            "npu" : "npu",
            "njets" : "NCleanedJets",
            "iso_1" : "iso_1",
            "iso_2" : "iso_2",
            "ptvis" : "Boson_Pt"
            }
            



    reader=csv.reader(open(config['inputFile'],"rb"),delimiter=',')
    datacsv=list(reader)
    header = np.array(datacsv[0]).astype(np.str)
    inputdatentot =np.array(datacsv[1:]).astype(np.float32)
    
    dictInputTot = {}
    for index in range(0,header.shape[0]):
        if header[index] in ArtusDict:
            dictInputTot[ArtusDict[header[index]]] = index
        else:
            dictInputTot[header[index]] = index
            #print(header[index])
    
    #for name in dictInputTot:
       	#print(name)


    plotnames = config[config['activePlotting']]['inputVariables']


    inputDataPlot = np.empty(shape=[inputdatentot.shape[0],0]).astype(np.float32)
    

    dictPlot = {}




    dt = int((time.time() - start))
    print('Elapsed time for loading dataset: ', dt)
    
    
    for index, entry in enumerate(dictInputTot):
	if entry in plotnames:
	    dictPlot[entry] = inputDataPlot.shape[1]
	    inputDataPlot = np.hstack((inputDataPlot, np.array(inputdatentot[:,dictInputTot[entry]]).reshape(inputdatentot.shape[0],1)))

    print(dictPlot)
    print(inputDataPlot.shape)


    
    dt = int((time.time() - start))
    print('Elapsed time for loading whole data: ', dt)

    return inputDataPlot, dictPlot
    
    
    


def main(config):
	    

	# Load the dataset
	print("Loading data...")
	plotData, dictPlot = load_datasetcsv(config)

	#constraints to data
	if 'constraints' in config:
	    print("Size of dataset before applying constraint: %i"%plotData.shape[0])
	    exec("{}".format("plotData = plotData[%s,:]"%config['constraints']))
	    print("Size of dataset after applying constraint: %i"%plotData.shape[0])

	#corrMatrix = np.corrcoef(plotData,rowvar=0)

	mean = plotData.mean()
	std = plotData.std()
	plotData = (plotData - mean)/std

	#print(corrMatrix.shape)
	results = PCA(plotData)

	print(results.fracs)

	print(results.Y[0])	
	print(results.Y[1])	
	print(results.Y[2])	
	print(np.sum(results.Y[0]*results.Y[1]))
	x = []
	y = []
	z = []
	for item in results.Y:
	    x.append(item[0])
	    y.append(item[1])
	    z.append(item[2])

	plt.close('all') # close all latent plotting windows
	fig1 = plt.figure() # Make a plotting figure
	ax = Axes3D(fig1) # use the plotting figure to create a Axis3D object.
	pltData = [x,y,z] 
	ax.scatter(pltData[0], pltData[1], pltData[2], 'bo') # make a scatter plot of blue dots from the data
	 
	# make simple, bare axis lines through space:
	xAxisLine = ((min(pltData[0]), max(pltData[0])), (0, 0), (0,0)) # 2 points make the x-axis line at the data extrema along x-axis 
	ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r') # make a red line for the x-axis.
	yAxisLine = ((0, 0), (min(pltData[1]), max(pltData[1])), (0,0)) # 2 points make the y-axis line at the data extrema along y-axis
	ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r') # make a red line for the y-axis.
	zAxisLine = ((0, 0), (0,0), (min(pltData[2]), max(pltData[2]))) # 2 points make the z-axis line at the data extrema along z-axis
	ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r') # make a red line for the z-axis.
	 
	# label the axes 
	ax.set_xlabel("x-axis label") 
	ax.set_ylabel("y-axis label")
	ax.set_zlabel("y-axis label")
	ax.set_title("The title of the plot")
	plt.savefig('correlation.png') # show the plot

	if 'export' in config:
	    bashCommand = 'cp -r %s /ekpwww/nzaeh/public_html/'%config['outputDir']
	    os.system(bashCommand)
	    print('Plots exported to ekpwww!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make MVAMet control plots.')
    parser.add_argument('-p', '--plottingconfig', default='../configs/config.json', help='Path to configurations file')
    parser.add_argument('-i', '--inputfile', default='', help='[optional] Inputfile from which to create the plots from')
    parser.add_argument('-o', '--outputfolder', default='', help='[optional] Foldername in which to store the plots in')
    parser.add_argument('-c', '--constraints', default='', help='[optional] Constraints to data. E.g.: 50<=plotData[:,dictPlot["Boson_Pt"]] & 50<=plotData[:,dictPlot["recoilslimmedMETs_LongZ"]]')
    parser.add_argument('-e', '--export', dest='export', action='store_true', help='[optional] Exports plots to ekpwww after creating them')
    parser.set_defaults(export=False)
    args = parser.parse_args()
    print('Used Configfile: ',args.plottingconfig)



    with open(args.plottingconfig, 'r') as f:
	config = json.load(f)
    
    if args.export:
	print('Exporting files to ekpwww afterwards.')
	config['export'] = True

    if not args.inputfile == '':
        config['inputFile'] = args.inputfile
    
    if not args.outputfolder == '':
        config['outputDir'] = "../plots/" + args.outputfolder + "/"
    
    if not args.constraints == '':
        config['constraints'] = args.constraints
    	print('Inputfile: %s with constraint %s'%(config['inputFile'],args.constraints))
    else:
    	print('Inputfile: %s'%config['inputFile'])
	
    
    main(config)


