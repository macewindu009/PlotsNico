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
from matplotlib.colors import LogNorm

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
            


    filename = config['inputFile'][0]
    print('Loading dataset %s ...'%filename)
    reader=csv.reader(open(filename,"rb"),delimiter=',')
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


    plotnames = config[config['activePlotting']]['plotVariables']


    inputDataPlot = np.empty(shape=[inputdatentot.shape[0],0]).astype(np.float32)
    

    dictPlot = {}




    dt = int((time.time() - start))
    print('Elapsed time for loading dataset: ', dt)
    lastTime = time.time()
    
    
    for index, entry in enumerate(dictInputTot):
	if entry in plotnames:
	    dictPlot[entry] = inputDataPlot.shape[1]
	    inputDataPlot = np.hstack((inputDataPlot, np.array(inputdatentot[:,dictInputTot[entry]]).reshape(inputdatentot.shape[0],1)))

    trainingheader =  ["LongZCorrectedRecoil_LongZ","LongZCorrectedRecoil_PerpZ","LongZCorrectedRecoil_Phi"]
    for index in range(len(config['inputFile'])-1):
    	filename = config['inputFile'][index+1]
	print('Loading dataset %s ...'%filename)
    	reader=csv.reader(open(filename,"rb"),delimiter=',')
    	datacsv=list(reader)
    	header = np.array(datacsv[0]).astype(np.str)
    	inputdatentot =np.array(datacsv[1:]).astype(np.float32)
    	for indexHeader in range(0,header.shape[0]):
            if header[indexHeader] in trainingheader:
		dictPlot['V' + str(index)+header[indexHeader]] = inputDataPlot.shape[1] 
		inputDataPlot = np.hstack((inputDataPlot, np.array(inputdatentot[:,indexHeader]).reshape(inputdatentot.shape[0],1)))
	dt = int(time.time() - lastTime)
        lastTime = time.time()
    	print('Elapsed time for loading dataset: ', dt)

    print(dictPlot)
    print(inputDataPlot.shape)

    
    dt = int((time.time() - start))
    print('Elapsed time for loading whole data: ', dt)

    return inputDataPlot, dictPlot
    
    
    
def make_ResponseCorrectedPlot(config, XRange, YStd, YResponse, bosonName, targetvariable,  minrange,maxrange, stepwidth, ptmin,ptmax, labelname = 'MVAMet', relateVar = 'p_t^Z', relateUnits = 'GeV'):

    plt.clf()
    ResCorr = YStd[:]/YResponse[:]
    plt.figure(6)
    plt.plot(XRange[:-1]+stepwidth/2.,ResCorr,'o-',label=labelname)
    plt.figure(0)
    plt.clf()    
    
    
    return 
    
def make_ResolutionPlot(config,plotData,dictPlot, bosonName, targetvariable,  minrange=42,maxrange=0, stepwidth=0, ptmin =0,ptmax=0, labelname = 'MVAMet', relateVar = 'p_t^Z', relateUnits = 'GeV'):


    if minrange == 42:
	minrange = plotData[:,targetindex].min()
    if maxrange == 0:
	maxrange = plotData[:,targetindex].max()
    if stepwidth == 0:
	stepwidth = (maxrange-minrange)/20
    XRange = np.arange(minrange,maxrange,stepwidth)
    YMean = np.zeros((XRange.shape[0]-1,1))
    YStd = np.zeros((XRange.shape[0]-1,1))
    
    print('Resolution %s versus %s'%(bosonName,targetvariable))
    #YValues 
    for index in range(0,XRange.shape[0]-1):

	AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]
	currentDistri = AlternativeDistri[:,dictPlot[bosonName]]+AlternativeDistri[:,dictPlot['Boson_Pt']]

	if currentDistri.shape == (0,1):
	    YMean[index] = 0
	    YStd[index] = 0
	else:
	    num_bins = 100
	    
            fitDistri = currentDistri[((currentDistri.mean()-4*currentDistri.std())<currentDistri[:]) & (currentDistri[:] <(currentDistri.mean()+4*currentDistri.std()))]
            
            YMean[index] = fitDistri.mean()
            YStd[index] = fitDistri.std()
                
        
                
    plt.figure(5)
    plt.plot(XRange[:-1]+stepwidth/2.,YStd[:],'o-',label=labelname)
    plt.figure(0)
    plt.clf()
    
    return XRange, YStd

   
def make_METResolutionPlot(config,plotData,dictPlot, bosonName, targetvariable,  minrange=42,maxrange=0, stepwidth=0, ptmin =0,ptmax=0, labelname = 'MVAMet', relateVar = 'p_t^Z', relateUnits = 'GeV'):


    #XRange = np.arange(plotData[:,targetindex].min(),plotData[:,targetindex].max(),(plotData[:,targetindex].max()-plotData[:,targetindex].min())/nbins)
    if minrange == 42:
	minrange = plotData[:,targetindex].min()
    if maxrange == 0:
	maxrange = plotData[:,targetindex].max()
    if stepwidth == 0:
	stepwidth = (maxrange-minrange)/20
    XRange = np.arange(minrange,maxrange,stepwidth)
    YMean = np.zeros((XRange.shape[0]-1,1))
    YStd = np.zeros((XRange.shape[0]-1,1))
    
    print('MET Resolution %s versus %s'%(bosonName,targetvariable))
    #YValues 
    for index in range(0,XRange.shape[0]-1):

	AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]
	currentDistri = AlternativeDistri[:,dictPlot[bosonName]]-AlternativeDistri[:,dictPlot['genMet_Pt']]

	if currentDistri.shape == (0,1):
	    YMean[index] = 0
	    YStd[index] = 0
	else:
	    fitDistri = currentDistri[((currentDistri.mean()-4*currentDistri.std())<currentDistri[:]) & (currentDistri[:] <(currentDistri.mean()+4*currentDistri.std()))]

            YMean[index] = fitDistri.mean()
            YStd[index] = fitDistri.std()
    plt.figure(10)
    plt.plot(XRange[:-1]+stepwidth/2.,YStd[:],'o-',label=labelname)
    plt.figure(0)
    plt.clf()
    
    
    return 0
    


def make_ResponsePlot(config, plotData,dictPlot, bosonName, targetvariable,  minrange=42,maxrange=0, stepwidth=0, ptmin=0, ptmax=0, labelname = 'MVAMet', relateVar = 'p_t^Z', relateUnits = 'GeV'):
  
    #XRange = np.arange(plotData[:,targetindex].min(),plotData[:,targetindex].max(),(plotData[:,targetindex].max()-plotData[:,targetindex].min())/nbins)
    if minrange == 42:
	minrange = plotData[:,dictPlot[targetvariable]].min()
    if maxrange == 0:
	maxrange = plotData[:,dictPlot[targetvariable]].max()
    if stepwidth == 0:
	stepwidth = (maxrange-minrange)/20
	
    XRange = np.arange(minrange,maxrange,stepwidth)

    YMean = np.zeros((XRange.shape[0]-1,1))
    YStd = np.zeros((XRange.shape[0]-1,1))
    print('Response %s versus %s'%(bosonName,targetvariable))
    
    
    
    #YValues 
    ignoredEntries = 0
    for index in range(0,XRange.shape[0]-1):

	AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]

	currentDistri = -AlternativeDistri[:,dictPlot[bosonName]]/AlternativeDistri[:,dictPlot['Boson_Pt']]
        fitDistri = currentDistri[((currentDistri.mean()-4*currentDistri.std())<currentDistri[:]) & (currentDistri[:] <(currentDistri.mean()+4*currentDistri.std()))]
	YMean[index] = fitDistri.mean()
	YStd[index] = fitDistri.std()
    plt.clf()
    plt.figure(4)
    plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o-',label=labelname)
    
    plt.figure(0)
    
    return YMean

    
def make_METResponsePlot(config, plotData,dictPlot, bosonName, targetvariable,  minrange=42,maxrange=0, stepwidth=0, ptmin=0, ptmax=0, labelname = 'MVAMet', relateVar = 'p_t^Z', relateUnits = 'GeV'):
  
    #XRange = np.arange(plotData[:,targetindex].min(),plotData[:,targetindex].max(),(plotData[:,targetindex].max()-plotData[:,targetindex].min())/nbins)
    if minrange == 42:
	minrange = plotData[:,dictPlot[targetvariable]].min()
    if maxrange == 0:
	maxrange = plotData[:,dictPlot[targetvariable]].max()
    if stepwidth == 0:
	stepwidth = (maxrange-minrange)/20
	
    XRange = np.arange(minrange,maxrange,stepwidth)
    YMean = np.zeros((XRange.shape[0]-1,1))
    YStd = np.zeros((XRange.shape[0]-1,1))
    print('MET Response %s versus %s'%(bosonName,targetvariable))
    
    #YValues 
    ignoredEntries = 0
    for index in range(0,XRange.shape[0]-1):

	AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]

	currentDistri = -AlternativeDistri[:,dictPlot[bosonName]]/AlternativeDistri[:,dictPlot['genMet_Pt']]
        fitDistri = currentDistri[((currentDistri.mean()-4*currentDistri.std())<currentDistri[:]) & (currentDistri[:] <(currentDistri.mean()+4*currentDistri.std()))]

        YMean[index] = fitDistri.mean()
        YStd[index] = fitDistri.std()
    plt.clf()
    plt.figure(9)
    plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o-',label=labelname)
    
    plt.figure(0)
    

    return 0

    
    
def make_ResolutionPerpPlot(config,plotData,dictPlot, bosonName, targetvariable,  minrange=42,maxrange=0, stepwidth=0, ptmin =0,ptmax=0, labelname = 'MVAMet', relateVar = 'p_t^Z', relateUnits = 'GeV'):


    if minrange == 42:
	minrange = plotData[:,targetindex].min()
    if maxrange == 0:
	maxrange = plotData[:,targetindex].max()
    if stepwidth == 0:
	stepwidth = (maxrange-minrange)/20
    XRange = np.arange(minrange,maxrange,stepwidth)
    YMean = np.zeros((XRange.shape[0]-1,1))
    YStd = np.zeros((XRange.shape[0]-1,1))
    
    print('Resolution Perp %s versus %s'%(bosonName,targetvariable))
    #YValues 
    for index in range(0,XRange.shape[0]-1):

	AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]
	currentDistri = AlternativeDistri[:,dictPlot[bosonName]]

	if currentDistri.shape == (0,1):
	    YMean[index] = 0
	    YStd[index] = 0
	else:
            fitDistri = currentDistri[((currentDistri.mean()-4*currentDistri.std())<currentDistri[:]) & (currentDistri[:] <(currentDistri.mean()+4*currentDistri.std()))]

            YMean[index] = fitDistri.mean()
            YStd[index] = fitDistri.std()
	
            
	

    plt.figure(8)
    plt.plot(XRange[:-1]+stepwidth/2.,YStd[:],'o-',label=labelname)
    plt.figure(0)
    plt.clf()
    return 
    
def make_ResponsePerpPlot(config, plotData,dictPlot, bosonName, targetvariable,  minrange=42,maxrange=0, stepwidth=0, ptmin=0, ptmax=0, labelname = 'MVAMet', relateVar = 'p_t^Z', relateUnits = 'GeV'):
  
    if minrange == 42:
	minrange = plotData[:,dictPlot[targetvariable]].min()
    if maxrange == 0:
	maxrange = plotData[:,dictPlot[targetvariable]].max()
    if stepwidth == 0:
	stepwidth = (maxrange-minrange)/20
	
    XRange = np.arange(minrange,maxrange,stepwidth)
    YMean = np.zeros((XRange.shape[0]-1,1))
    YStd = np.zeros((XRange.shape[0]-1,1))
    print('Response Perp %s versus %s'%(bosonName,targetvariable))
    
    
    #YValues 
    ignoredEntries = 0
    for index in range(0,XRange.shape[0]-1):

	AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]

	currentDistri = AlternativeDistri[:,dictPlot[bosonName]]
        fitDistri = currentDistri[((currentDistri.mean()-4*currentDistri.std())<currentDistri[:]) & (currentDistri[:] <(currentDistri.mean()+4*currentDistri.std()))]

        YMean[index] = fitDistri.mean()
        YStd[index] = fitDistri.std()
    plt.clf()
    plt.figure(7)
    plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o-',label=labelname)
    plt.figure(0)
    
    return 

    

def make_ControlPlots(config, plotData,dictPlot, bosonName, targetvariable, minrange=42,maxrange=0, stepwidth=0, ptmin=0,ptmax=0, labelname = 'MVAMet', relateVar = 'p_t^Z', relateUnits = 'GeV'):

    bosonNameLong = bosonName + '_LongZ'
    bosonNamePerp = bosonName + '_PerpZ'
    maxrange += stepwidth
    if not os.path.exists((config['outputDir'] + 'ControlPlots/SingleDistributions/')):
	os.makedirs((config['outputDir'] + 'ControlPlots/SingleDistributions/'))
    XRange, YVariance = make_ResolutionPlot(config, plotData, dictPlot, bosonNameLong, targetvariable,  minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits) 
    YResponse = make_ResponsePlot(config, plotData, dictPlot, bosonNameLong, targetvariable,  minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits)
    make_ResponseCorrectedPlot(config, XRange, YVariance, YResponse, bosonNameLong, targetvariable,  minrange,maxrange, stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits)
    make_ResolutionPerpPlot(config, plotData, dictPlot, bosonNamePerp, targetvariable,  minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits) 
    make_ResponsePerpPlot(config, plotData, dictPlot, bosonNamePerp, targetvariable,  minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits)
    
    if bosonName == "LongZCorrectedRecoil" and "recoMetOnGenMetProjectionPar" in dictPlot and not plotData[dictPlot["recoMetOnGenMetProjectionPar"]].mean() == -999:
        make_METResponsePlot(config, plotData, dictPlot, "recoMetOnGenMetProjectionPar", targetvariable,  minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits)
        make_METResolutionPlot(config, plotData, dictPlot, "recoMetOnGenMetProjectionPar", targetvariable,  minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits) 
    
    if bosonName == "recoilslimmedMETs" and "recoPfMetOnGenMetProjectionPar" in dictPlot and not plotData[dictPlot["recoPfMetOnGenMetProjectionPar"]].mean() == -999:
        make_METResponsePlot(config, plotData, dictPlot, "recoPfMetOnGenMetProjectionPar", targetvariable,  minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits)
        make_METResolutionPlot(config, plotData, dictPlot, "recoPfMetOnGenMetProjectionPar", targetvariable,  minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits) 
    
    
    
    return 
    
    
def plot_results(config, plotData, dictPlot):
    
    #plotData = plotData[0==plotData[:,dictPlot['NCleanedJets']],:]
    
    
    
    
    
    
    
    
    plotconfig = config[config['activePlotting']]
    
    num_bins = 50
  
    if not os.path.exists(config['outputDir']):
	os.makedirs(config['outputDir'])
	
    
    bosonmin = [0,0,plotconfig['BosonCut']]
    bosonmax = [plotData[:,dictPlot['Boson_Pt']].max(),plotconfig['BosonCut'],plotData[:,dictPlot['Boson_Pt']].max()]
    
    #BDT Performance Plots 
#Phi ControlPlots
    for i, min in enumerate(bosonmin):
        slicedData = plotData[min<=plotData[:,dictPlot['Boson_Pt']],:]
        slicedData = slicedData[bosonmax[i]>=slicedData[:,dictPlot['Boson_Pt']],:]
        if plotconfig['plotBDTPerformance']:
            #BDT
            if 'LongZCorrectedRecoil_LongZ' in dictPlot:
		if 'inputLabel' in config:
                    make_ControlPlots(config, slicedData, dictPlot, 'LongZCorrectedRecoil', 'NVertex', 5,40,5,min,bosonmax[i],config['inputLabel'][0],'\# PV','')
		else:
                    make_ControlPlots(config, slicedData, dictPlot, 'LongZCorrectedRecoil', 'NVertex', 5,40,5,min,bosonmax[i],'MVAMet 1','\# PV','')
            
             
    	    for index in range(len(config['inputFile'])-1):
            	if 'V%iLongZCorrectedRecoil_LongZ'%index in dictPlot:
		    if 'inputLabel' in config:
                     	make_ControlPlots(config, slicedData, dictPlot, 'V%iLongZCorrectedRecoil'%index, 'NVertex', 5,40,5,min,bosonmax[i],config['inputLabel'][index+1],'\# PV','')
		    else:
                     	make_ControlPlots(config, slicedData, dictPlot, 'V%iLongZCorrectedRecoil'%index, 'NVertex', 5,40,5,min,bosonmax[i],'MVAMet %i'%(index+2),'\# PV','')
        
        
        #slimmedMet
        if 'recoilslimmedMETs_LongZ' in dictPlot:
            make_ControlPlots(config, slicedData, dictPlot, 'recoilslimmedMETs', 'NVertex', 5,40,5,min,bosonmax[i],'PFMet','\# PV','')
        
        
        
       
       
	plt.clf()
	plt.figure(4) 
        plt.xlabel(r'#PV ($%i\,\mathrm{GeV} < p_t^Z < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]), fontsize = 17)
    	plt.ylabel(r'$<U_\| / p_t^Z>$',fontsize = 17)
    	plt.title('Response 'r'$U_\|$'' vs 'r'#PV',fontsize = 20)
        legend = plt.legend(loc='best', shadow=True)
    	plt.plot((plt.xlim()[0], plt.xlim()[1]), (1, 1), 'k--')
        plt.savefig(config['outputDir'] + 'Response_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()

        plt.figure(5)
        plt.xlabel(r'#PV ($%i\,\mathrm{GeV} < p_t^Z < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]),fontsize = 17)
    	plt.ylabel(r'$\sigma(<U_\| - p_t^Z>)$'' in GeV',fontsize = 17)
    	plt.title('Resolution 'r'$U_\|$'' vs 'r'#PV',fontsize = 20)
        legend = plt.legend(loc='best', shadow=True)
        plt.savefig(config['outputDir'] + 'Resolution_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()

        plt.figure(6)
        plt.xlabel(r'#PV ($%i\,\mathrm{GeV} < p_t^Z < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]), fontsize = 17)
        plt.ylabel('Resolution / Response',fontsize = 17)
    	plt.title('Response Corrected vs 'r'#PV',fontsize = 20)
        legend = plt.legend(loc='best', shadow=True)
        plt.savefig(config['outputDir'] + 'ResponseCorrected_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()

        plt.figure(7)
        legend = plt.legend(loc='best', shadow=True)
        plt.xlabel(r'#PV ($%i\,\mathrm{GeV} < p_t^Z < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]), fontsize = 17)
    	plt.ylabel(r'$<U_\bot>$',fontsize = 17)
    	plt.title('Response 'r'$U_\bot$'' vs 'r'#PV',fontsize = 20)
        plt.savefig(config['outputDir'] + 'ResponsePerp_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()

        plt.figure(8)
        plt.xlabel(r'#PV ($%i\,\mathrm{GeV} < p_t^Z < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]), fontsize = 17)
    	plt.ylabel(r'$\sigma(<U_\bot>)$',fontsize = 17)
    	plt.title('Resolution 'r'$U_\bot$'' vs 'r'#PV',fontsize = 20)
        legend = plt.legend(loc='best', shadow=True)
        plt.savefig(config['outputDir'] + 'ResolutionPerp_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()

        plt.figure(9)
        plt.xlabel(r'#PV ($%i\,\mathrm{GeV} < p_t^Z < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]), fontsize = 17)
        plt.ylabel(r'$E_{t\|}^{miss}/E_{t,gen}^{miss}$',fontsize = 17)
        legend = plt.legend(loc='best', shadow=True)
        plt.savefig(config['outputDir'] + 'METResponse_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()

        plt.figure(10)
        plt.ylabel('std(MET_Long/genMET)',fontsize = 17)
        plt.ylabel(r'$\sigma(E_{t\|}^{miss}-E_{t,gen}^{miss})$',fontsize = 17)
        plt.xlabel(r'#PV ($%i\,\mathrm{GeV} < p_t^Z < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]), fontsize = 17)
        legend = plt.legend(loc='best', shadow=True)
        plt.savefig(config['outputDir'] + 'METResolution_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()
        plt.figure(0)
        

    
    #Boson PT
    if plotconfig['plotBDTPerformance']:
        if 'LongZCorrectedRecoil_LongZ' in dictPlot:
	    if 'inputLabel' in config:
            	make_ControlPlots(config, plotData, dictPlot, 'LongZCorrectedRecoil', 'Boson_Pt', 10,200,10,0,0,config['inputLabel'][0],'p_t^Z','GeV')
	    else:
            	make_ControlPlots(config, plotData, dictPlot, 'LongZCorrectedRecoil', 'Boson_Pt', 10,200,10,0,0,'MVAMet 1','p_t^Z','GeV')

    	for index in range(len(config['inputFile'])-1):
    	    if 'V%iLongZCorrectedRecoil_LongZ'%index in dictPlot:
		if 'inputLabel' in config:
        	    make_ControlPlots(config, plotData, dictPlot, 'V%iLongZCorrectedRecoil'%index, 'Boson_Pt', 10,200,10,0,0,config['inputLabel'][index+1],'p_t^Z','GeV')
		else:
        	    make_ControlPlots(config, plotData, dictPlot, 'V%iLongZCorrectedRecoil'%index, 'Boson_Pt', 10,200,10,0,0,'MVAMet %i'%(index+2),'p_t^Z','GeV')

    
    

    if 'recoilslimmedMETs_LongZ' in dictPlot:
        make_ControlPlots(config, plotData, dictPlot, 'recoilslimmedMETs', 'Boson_Pt', 10,200,10,0,0,'PFMet','p_t^Z','GeV')
    
    
    


    plt.figure(4)
    legend = plt.legend(loc='best', shadow=True)
    #legend.get_frame().set_alpha(0.5)
    plt.xlabel(r'$p_t^Z$'' in GeV', fontsize = 17)
    plt.ylabel(r'$<U_\| / p_t^Z>$',fontsize = 17)
    plt.title('Response 'r'$U_\|$'' vs 'r'$p_t^Z$',fontsize = 20)
    plt.plot((plt.xlim()[0], plt.xlim()[1]), (1, 1), 'k--')
    plt.savefig(config['outputDir'] + 'Response_vs_BosonPt')
    plt.clf()
    plt.figure(5)
    plt.xlabel(r'$p_t^Z$'' in GeV', fontsize = 17)
    plt.ylabel(r'$\sigma(<U_\| - p_t^Z>)$'' in GeV',fontsize = 17)
    plt.title('Resolution 'r'$U_\|$'' vs 'r'$p_t^Z$',fontsize = 20)
    legend = plt.legend(loc='best', shadow=True)
    #legend.get_frame().set_alpha(0.5)
    plt.savefig(config['outputDir'] + 'Resolution_vs_BosonPt')
    plt.clf()
    plt.figure(6)
    plt.xlabel(r'$p_t^Z$'' in GeV', fontsize = 17)
    plt.ylabel('Resolution / Response',fontsize = 17)
    plt.title('Response Corrected vs 'r'$p_t^Z$',fontsize = 20)
    legend = plt.legend(loc='best', shadow=True)
    #legend.get_frame().set_alpha(0.5)
    plt.savefig(config['outputDir'] + 'ResponseCorrected_vs_BosonPt')
    plt.clf()
    plt.figure(7)
    legend = plt.legend(loc='best', shadow=True)
    #legend.get_frame().set_alpha(0.5)
    plt.xlabel(r'$p_t^Z$'' in GeV', fontsize = 17)
    plt.ylabel(r'$<U_\bot>$',fontsize = 17)
    plt.title('Response 'r'$U_\bot$'' vs 'r'$p_t^Z$',fontsize = 20)
    plt.savefig(config['outputDir'] + 'ResponsePerp_vs_BosonPt')
    plt.clf()
    plt.figure(8)
    plt.xlabel(r'$p_t^Z$'' in GeV', fontsize = 17)
    plt.ylabel(r'$\sigma(<U_\bot>)$',fontsize = 17)
    plt.title('Resolution 'r'$U_\bot$'' vs 'r'$p_t^Z$',fontsize = 20)
    legend = plt.legend(loc='best', shadow=True)
    #legend.get_frame().set_alpha(0.5)
    plt.savefig(config['outputDir'] + 'ResolutionPerp_vs_BosonPt')
    plt.clf()
    plt.figure(9)
    legend = plt.legend(loc='best', shadow=True)
    #legend.get_frame().set_alpha(0.5)
    plt.xlabel(r'$p_t^Z$'' in GeV', fontsize = 17)
    #plt.ylabel('MET_Long/genMET',fontsize = 17)
    plt.ylabel(r'$E_{t\|}^{miss}/E_{t,gen}^{miss}$',fontsize = 17)
    #plt.ylabel(r'$\ensuremath{{\not\mathrel{E}}_T}$',fontsize = 17)
    plt.savefig(config['outputDir'] + 'METResponse_vs_BosonPt')
    plt.clf()
    plt.figure(10)
    plt.xlabel(r'$p_t^Z$'' in GeV', fontsize = 17)
    plt.ylabel(r'$\sigma(E_{t\|}^{miss}-E_{t,gen}^{miss})$',fontsize = 17)
    legend = plt.legend(loc='best', shadow=True)
    #legend.get_frame().set_alpha(0.5)
    plt.savefig(config['outputDir'] + 'METResolution_vs_BosonPt')
    plt.clf()
    plt.figure(0)
	    
    

    print('Plots created')
    
    return True

    



def main(config):
    

    # Load the dataset
    print("Loading data...")
    plotData, dictPlot = load_datasetcsv(config)
    

    #constraints to data
    if 'constraints' in config:
        print("Size of dataset before applying constraint: %i"%plotData.shape[0])
        exec("{}".format("plotData = plotData[%s,:]"%config['constraints']))
        print("Size of dataset after applying constraint: %i"%plotData.shape[0])



    print("plotDatashape", plotData.shape)
    plot_results(config, plotData, dictPlot)

    if 'export' in config:
        bashCommand = 'cp -r %s /ekpwww/nzaeh/public_html/'%config['outputDir']
        os.system(bashCommand)
        print('Plots exported to ekpwww!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make MVAMet control plots.')
    parser.add_argument('-p', '--plottingconfig', default='../configs/config.json', help='Path to configurations file')
    parser.add_argument('-i', '--inputfile',nargs='+', default='', help='[optional] Inputfile from which to create the plots from')
    parser.add_argument('-n', '--inputlabels',nargs='+', default='', help='[optional] Inputlabelname to use in plots')
    parser.add_argument('-o', '--outputfolder', default='', help='[optional] Foldername in which to store the plots in')
    parser.add_argument('-c', '--constraints', default='', help='[optional] Constraints to data. E.g.: 50<=plotData[:,dictPlot["Boson_Pt"]] & 50<=plotData[:,dictPlot["recoilslimmedMETs_LongZ"]]')
    parser.add_argument('-e', '--export', dest='export', action='store_true', help='[optional] Exports plots to ekpwww after creating them')
    parser.set_defaults(export=False)
    args = parser.parse_args()
    print('Used Configfile: ',args.plottingconfig)



    with open(args.plottingconfig, 'r') as f:
	config = json.load(f)
    
    if not args.inputfile == '':
        config['inputFile'] = args.inputfile
    	
    if not args.inputfile == '':
        config['inputLabel'] = args.inputlabels

    print('Inputfile: ')
    for names in config['inputFile']:
	print(names)	

    if not args.constraints == '':
        config['constraints'] = args.constraints
	print('Applying constraints %s'%args.constraints)

    if args.export:
	print('Plots will be exported to ekpwww afterwards.') 
	config['export'] = True

    if not args.outputfolder == '':
        config['outputDir'] = "../plots/" + args.outputfolder + "/"
    
    
    
    main(config)


