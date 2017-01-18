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
    
    for name in dictInputTot:
        print(name)


    plotnames = config[config['activePlotting']]['plotVariables']


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
    
    
    
def make_Plot(variablename, inputData, dictPlot, outputdir):

    histData = inputData[:,dictPlot[variablename]]
    
    if not os.path.exists(outputdir):
	os.makedirs(outputdir)
	
    
    
    num_bins = 100
    
    if variablename == 'targetRecoilFromBDT' or variablename == 'targetRecoilFromSlimmed':
	n, bins, patches = plt.hist(histData, num_bins, facecolor='green', alpha=0.5, range=[-50, 50])
    else:
	n, bins, patches = plt.hist(histData, num_bins, facecolor='green', alpha=0.5)
    plt.xlabel(variablename)
    plt.ylabel('Hits')

	
    plt.savefig((outputdir+variablename+".png"))
    plt.clf()
    return 0
    
    
def make_ResponseCorrectedPlot(config, XRange, YStd, YResponse, bosonName, targetvariable, resultData, dictResult, minrange,maxrange, stepwidth, ptmin,ptmax):

    plt.clf()
    ResCorr = YStd[:]/YResponse[:]
    plt.plot(XRange[:-1]+stepwidth/2.,ResCorr,'o')
    plt.xlabel(targetvariable)
    plt.ylabel('Resolution / Response')
    if ptmax == 0:
	  plt.savefig(config['outputDir'] + "ControlPlots/ResponseCorrected_%s_vs_%s.png" %(bosonName,targetvariable))
    else:
	  plt.savefig(config['outputDir'] + "ControlPlots/ResponseCorrected(%i<Pt<%i)_%s_vs_%s.png" %(ptmin,ptmax,bosonName,targetvariable))
    plt.figure(6)
    plt.plot(XRange[:-1]+stepwidth/2.,ResCorr,'o',label=bosonName)
    plt.figure(0)
    plt.clf()    
    
    
    if bosonName == 'LongZCorrectedRecoil_LongZ':
    #if bosonName == 'recoilslimmedMETs_LongZ':
	if ptmax > config[config['activePlotting']]['BosonCut']:
	    Upper = 'Max'
	else:
	    Upper = 'Cut'
	    
	if ptmin < config[config['activePlotting']]['BosonCut']:
	    Lower = 'Min'
	else:
	    Lower = 'Cut'
      
	if targetvariable == 'Boson_Pt':
	    resultvalue = ResCorr[np.isfinite(ResCorr)].min()
	    dictResult["ResCorr_Min_%s_MintoMax"%(targetvariable)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
	    XRangeTemp = XRange[:-1]
	    
	    resultvalue = ResCorr[np.isfinite(ResCorr)].sum()/ResCorr[np.isfinite(ResCorr)].shape[0]
	    dictResult["ResCorr_Int_%s_MintoMax"%(targetvariable)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
	    XRangeTemp = XRange[:-1]
	    
	    resultvalue = ResCorr[np.isfinite(ResCorr)].max()
	    dictResult["ResCorr_Max_%s_MintoMax"%(targetvariable)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
	    XRangeTemp = XRange[:-1]
	    
	elif targetvariable == 'NVertex':

	    resultvalue = ResCorr[np.isfinite(ResCorr)].min()
	    dictResult["ResCorr_Min_%s_%sto%s"%(targetvariable,Lower,Upper)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
	    
	    resultvalue = ResCorr[np.isfinite(ResCorr)].sum()/ResCorr[np.isfinite(ResCorr)].shape[0]
	    dictResult["ResCorr_Int_%s_%sto%s"%(targetvariable,Lower,Upper)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
	
	    resultvalue = ResCorr[np.isfinite(ResCorr)].max()
	    dictResult["ResCorr_Max_%s_%sto%s"%(targetvariable,Lower,Upper)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
    
    
    
    return resultData, dictResult
    
def make_ResolutionPlot(config,plotData,dictPlot, bosonName, targetvariable, resultData, dictResult, minrange=42,maxrange=0, stepwidth=0, ptmin =0,ptmax=0):


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
    
    print('Resolution %s versus %s'%(bosonName,targetvariable))
    #YValues 
    for index in range(0,XRange.shape[0]-1):

	AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]
	currentDistri = AlternativeDistri[:,dictPlot[bosonName]]+AlternativeDistri[:,dictPlot['Boson_Pt']]

	if currentDistri.shape == (0,1):
	    YMean[index] = 0
	    YStd[index] = 0
	else:
	    """
	    print('current bin ',index,' entries:',currentDistri.shape[0])	
	    print('current bin ',index,' min:',currentDistri.min())
	    print('current bin ',index,' max:',currentDistri.max())
	    print('current bin ',index,' mean:',currentDistri.mean())
	    print('current bin ',index,' std:',currentDistri.std())
	    """
	    YMean[index] = currentDistri.mean()
	    YStd[index] = currentDistri.std()
	
	    if index < 12:
		plt.clf()
		num_bins = 50
		n, bins, patches = plt.hist(currentDistri, num_bins, normed=1, facecolor='green', alpha=0.5)
		y = mlab.normpdf(bins, currentDistri.mean(), currentDistri.std())
		plt.xlabel('%s at %f, mean: %f'%(targetvariable,(XRange[index+1]+XRange[index])/2,currentDistri.mean()))
		plt.ylabel('(MET Boson PT_Long) - (True Boson Pt)')
		plt.plot(bins, y, 'r--')
		if ptmax == 0:
		    plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/Resolution_%s_vs_%s_%i.png' %(bosonName,targetvariable, index)))
		else:
		    plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/Resolution(%i<Pt<%i)_%s_vs_%s_%i.png' %(ptmin,ptmax,bosonName,targetvariable, index)))
		    
    plt.clf()
    plt.plot(XRange[:-1]+stepwidth/2.,YStd[:],'o')
    plt.ylabel('(MET Boson PT_Long) - (True Boson Pt)')
    plt.xlabel(targetvariable)
    if ptmax == 0:
	plt.savefig(config['outputDir'] + "ControlPlots/Resolution_%s_vs_%s.png" %(bosonName,targetvariable))
    else:
	plt.savefig(config['outputDir'] + "ControlPlots/Resolution(%i<Pt<%i)_%s_vs_%s.png" %(ptmin,ptmax,bosonName,targetvariable))
    plt.figure(5)
    plt.plot(XRange[:-1]+stepwidth/2.,YStd[:],'o',label=bosonName)
    plt.figure(0)
    plt.clf()
    
    
    if bosonName == 'LongZCorrectedRecoil_LongZ':
    #if bosonName == 'recoilslimmedMETs_LongZ':
	if ptmax > config[config['activePlotting']]['BosonCut']:
	    Upper = 'Max'
	else:
	    Upper = 'Cut'
	    
	if ptmin < config[config['activePlotting']]['BosonCut']:
	    Lower = 'Min'
	else:
	    Lower = 'Cut'
      
	if targetvariable == 'Boson_Pt':
	    resultvalue = YStd[np.isfinite(YStd)].min()
	    dictResult["Resolution_Min_%s_MintoMax"%(targetvariable)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
	    XRangeTemp = XRange[:-1]
	    
	    resultvalue = YStd[np.isfinite(YStd)].sum()/YStd[np.isfinite(YStd)].shape[0]
	    dictResult["Resolution_Int_%s_MintoMax"%(targetvariable)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
	    XRangeTemp = XRange[:-1]
	    
	elif targetvariable == 'NVertex':

	    resultvalue = YStd[np.isfinite(YStd)].min()
	    dictResult["Resolution_Min_%s_%sto%s"%(targetvariable,Lower,Upper)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
	    
	    resultvalue = YStd[np.isfinite(YStd)].sum()/YStd[np.isfinite(YStd)].shape[0]
	    dictResult["Resolution_Int_%s_%sto%s"%(targetvariable,Lower,Upper)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))

    
    
    
    return resultData, dictResult, XRange, YStd

   
def make_METResolutionPlot(config,plotData,dictPlot, bosonName, targetvariable, resultData, dictResult, minrange=42,maxrange=0, stepwidth=0, ptmin =0,ptmax=0):


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
	currentDistri = AlternativeDistri[:,dictPlot[bosonName]]+AlternativeDistri[:,dictPlot['genMet_Pt']]

	if currentDistri.shape == (0,1):
	    YMean[index] = 0
	    YStd[index] = 0
	else:
	    """
	    print('current bin ',index,' entries:',currentDistri.shape[0])	
	    print('current bin ',index,' min:',currentDistri.min())
	    print('current bin ',index,' max:',currentDistri.max())
	    print('current bin ',index,' mean:',currentDistri.mean())
	    print('current bin ',index,' std:',currentDistri.std())
	    """
	    YMean[index] = currentDistri.mean()
	    YStd[index] = currentDistri.std()
	
	    if index < 12:
		plt.clf()
		num_bins = 50
		n, bins, patches = plt.hist(currentDistri, num_bins, normed=1, facecolor='green', alpha=0.5)
		y = mlab.normpdf(bins, currentDistri.mean(), currentDistri.std())
		plt.xlabel('%s at %f, mean: %f'%(targetvariable,(XRange[index+1]+XRange[index])/2,currentDistri.mean()))
		plt.ylabel('(MET) - (gen MET)')
		plt.plot(bins, y, 'r--')
		if ptmax == 0:
		    plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/METResolution_%s_vs_%s_%i.png' %(bosonName,targetvariable, index)))
		else:
		    plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/METResolution(%i<Pt<%i)_%s_vs_%s_%i.png' %(ptmin,ptmax,bosonName,targetvariable, index)))
		    
    plt.clf()
    plt.plot(XRange[:-1]+stepwidth/2.,YStd[:],'o')
    plt.ylabel('(MET) - (gen MET)')
    plt.xlabel(targetvariable)
    if ptmax == 0:
	plt.savefig(config['outputDir'] + "ControlPlots/METResolution_%s_vs_%s.png" %(bosonName,targetvariable))
    else:
	plt.savefig(config['outputDir'] + "ControlPlots/METResolution(%i<Pt<%i)_%s_vs_%s.png" %(ptmin,ptmax,bosonName,targetvariable))
    plt.figure(10)
    plt.plot(XRange[:-1]+stepwidth/2.,YStd[:],'o',label=bosonName)
    plt.figure(0)
    plt.clf()
    
    
    return 0
    
def make_ResponsePlot(config, plotData,dictPlot, bosonName, targetvariable, resultData, dictResult, minrange=42,maxrange=0, stepwidth=0, ptmin=0, ptmax=0):
  
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
	"""
	if AlternativeDistri.shape[0] == 0:
	    YMean[index] = 0
	    YStd[index] = 0
	else:
	"""
	YMean[index] = currentDistri.mean()
	YStd[index] = currentDistri.std()
	if index < 12:
	    plt.clf()
            num_bins = 50
	    n, bins, patches = plt.hist(currentDistri, num_bins, normed=1, facecolor='green', alpha=0.5)
	    y = mlab.normpdf(bins, currentDistri.mean(), currentDistri.std())
	    plt.xlabel('%s at %f, mean: %f'%(targetvariable,(XRange[index+1]+XRange[index])/2,currentDistri.mean()))
	    plt.ylabel('(MET Boson PT_Long)/(True Boson Pt)')
	    plt.plot(bins, y, 'r--')
	    if ptmax == 0:
		plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/Response_%s_vs_%s_%i.png' %(bosonName,targetvariable, index)))
	    else:
		plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/Response(%i<Pt<%i)_%s_vs_%s_%i.png' %(ptmin,ptmax,bosonName,targetvariable, index)))
    plt.clf()
    plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o')

    plt.xlabel(targetvariable)
    if ptmax == 0:
	plt.savefig(config['outputDir'] + "ControlPlots/Response_%s_vs_%s.png" %(bosonName,targetvariable))
    else:
	plt.savefig(config['outputDir'] + "ControlPlots/Response(%i<Pt<%i)_%s_vs_%s.png" %(ptmin,ptmax,bosonName,targetvariable))
    plt.clf()
    plt.figure(4)
    plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o',label=bosonName)
    
    plt.figure(0)
    

    if bosonName == 'LongZCorrectedRecoil_LongZ':
    #if bosonName == 'recoilslimmedMETs_LongZ':
	if ptmax > config[config['activePlotting']]['BosonCut']:
	    Upper = 'Max'
	else:
	    Upper = 'Cut'
	    
	if ptmin < config[config['activePlotting']]['BosonCut']:
	    Lower = 'Min'
	else:
	    Lower = 'Cut'
      
      
	if targetvariable == 'Boson_Pt':
	    resultvalue = ((YMean[np.isfinite(YMean)]-1)**2).sum()/YMean[np.isfinite(YMean)].shape[0]
	    dictResult["Response_Chi_%s_MintoMax"%(targetvariable)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
	    XRangeTemp = XRange[:-1]
	    
	    
	    YMeanCut = YMean[config[config['activePlotting']]['BosonCut']<XRangeTemp[:]]
	    print('shape cut:', YMeanCut.shape)
	    print('shape whole:', YMean.shape)
	    print('shape finite:', YMean[np.isfinite(YMean)].shape)
	    resultvalue = ((YMeanCut[np.isfinite(YMeanCut)]-1)**2).sum()/YMeanCut[np.isfinite(YMeanCut)].shape[0]
	    dictResult["Response_Chi_%s_CuttoMax"%(targetvariable)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
	elif targetvariable == 'NVertex':
	    
	    resultvalue = ((YMean[np.isfinite(YMean)]-1)**2).sum()/YMean[np.isfinite(YMean)].shape[0]
	    dictResult["Response_Chi_%s_%sto%s"%(targetvariable,Lower,Upper)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
    
    return resultData, dictResult, YMean

    
def make_METResponsePlot(config, plotData,dictPlot, bosonName, targetvariable, resultData, dictResult, minrange=42,maxrange=0, stepwidth=0, ptmin=0, ptmax=0):
  
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
	"""
	if AlternativeDistri.shape[0] == 0:
	    YMean[index] = 0
	    YStd[index] = 0
	else:
	"""
	YMean[index] = currentDistri.mean()
	YStd[index] = currentDistri.std()
	if index < 12:
	    plt.clf()
            num_bins = 50
	    n, bins, patches = plt.hist(currentDistri, num_bins, normed=1, facecolor='green', alpha=0.5)
	    y = mlab.normpdf(bins, currentDistri.mean(), currentDistri.std())
	    plt.xlabel('%s at %f, mean: %f'%(targetvariable,(XRange[index+1]+XRange[index])/2,currentDistri.mean()))
	    plt.ylabel('(MET)/(gen MET)')
	    plt.plot(bins, y, 'r--')
	    if ptmax == 0:
		plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/METResponse_%s_vs_%s_%i.png' %(bosonName,targetvariable, index)))
	    else:
		plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/METResponse(%i<Pt<%i)_%s_vs_%s_%i.png' %(ptmin,ptmax,bosonName,targetvariable, index)))
    plt.clf()
    plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o')

    plt.xlabel(targetvariable)
    if ptmax == 0:
	plt.savefig(config['outputDir'] + "ControlPlots/METResponse_%s_vs_%s.png" %(bosonName,targetvariable))
    else:
	plt.savefig(config['outputDir'] + "ControlPlots/METResponse(%i<Pt<%i)_%s_vs_%s.png" %(ptmin,ptmax,bosonName,targetvariable))
    plt.clf()
    plt.figure(9)
    plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o',label=bosonName)
    
    plt.figure(0)
    

    return 0

    
    
def make_ResolutionPerpPlot(config,plotData,dictPlot, bosonName, targetvariable, resultData, dictResult, minrange=42,maxrange=0, stepwidth=0, ptmin =0,ptmax=0):


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
    
    print('Resolution Perp %s versus %s'%(bosonName,targetvariable))
    #YValues 
    for index in range(0,XRange.shape[0]-1):

	AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]
	currentDistri = AlternativeDistri[:,dictPlot[bosonName]]

	if currentDistri.shape == (0,1):
	    YMean[index] = 0
	    YStd[index] = 0
	else:
	    """
	    print('current bin ',index,' entries:',currentDistri.shape[0])	
	    print('current bin ',index,' min:',currentDistri.min())
	    print('current bin ',index,' max:',currentDistri.max())
	    print('current bin ',index,' mean:',currentDistri.mean())
	    print('current bin ',index,' std:',currentDistri.std())
	    """
	    YMean[index] = currentDistri.mean()
	    YStd[index] = currentDistri.std()
	
	    if index < 12:
		plt.clf()
		num_bins = 50
		n, bins, patches = plt.hist(currentDistri, num_bins, normed=1, facecolor='green', alpha=0.5)
		y = mlab.normpdf(bins, currentDistri.mean(), currentDistri.std())
		plt.xlabel('%s at %f, mean: %f'%(targetvariable,(XRange[index+1]+XRange[index])/2,currentDistri.mean()))
		plt.ylabel('MET Boson PT_Perp')
		plt.plot(bins, y, 'r--')
		if ptmax == 0:
		    plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/ResolutionPerp_%s_vs_%s_%i.png' %(bosonName,targetvariable, index)))
		else:
		    plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/ResolutionPerp(%i<Pt<%i)_%s_vs_%s_%i.png' %(ptmin,ptmax,bosonName,targetvariable, index)))
		    
    plt.clf()
    plt.plot(XRange[:-1]+stepwidth/2.,YStd[:],'o')
    plt.ylabel('MET Boson PT_Perp')
    plt.xlabel(targetvariable)
    if ptmax == 0:
	plt.savefig(config['outputDir'] + "ControlPlots/ResolutionPerp_%s_vs_%s.png" %(bosonName,targetvariable))
    else:
	plt.savefig(config['outputDir'] + "ControlPlots/ResolutionPerp(%i<Pt<%i)_%s_vs_%s.png" %(ptmin,ptmax,bosonName,targetvariable))
    plt.figure(8)
    plt.plot(XRange[:-1]+stepwidth/2.,YStd[:],'o',label=bosonName)
    plt.figure(0)
    plt.clf()
    

    if bosonName == 'LongZCorrectedRecoil_PerpZ':
    #if bosonName == 'recoilslimmedMETs_LongZ':
	if ptmax > config[config['activePlotting']]['BosonCut']:
	    Upper = 'Max'
	else:
	    Upper = 'Cut'
	    
	if ptmin < config[config['activePlotting']]['BosonCut']:
	    Lower = 'Min'
	else:
	    Lower = 'Cut'
      
	if targetvariable == 'Boson_Pt':
	    resultvalue = YStd[np.isfinite(YStd)].min()
	    dictResult["ResolutionPerp_Min_%s_MintoMax"%(targetvariable)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
	    XRangeTemp = XRange[:-1]
	    
	    resultvalue = YStd[np.isfinite(YStd)].sum()/YStd[np.isfinite(YStd)].shape[0]
	    dictResult["ResolutionPerp_Int_%s_MintoMax"%(targetvariable)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
	    XRangeTemp = XRange[:-1]
	    
	elif targetvariable == 'NVertex':

	    resultvalue = YStd[np.isfinite(YStd)].min()
	    dictResult["ResolutionPerp_Min_%s_%sto%s"%(targetvariable,Lower,Upper)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
	    
	    resultvalue = YStd[np.isfinite(YStd)].sum()/YStd[np.isfinite(YStd)].shape[0]
	    dictResult["ResolutionPerp_Int_%s_%sto%s"%(targetvariable,Lower,Upper)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))

    
    
    
    return resultData, dictResult
    
def make_ResponsePerpPlot(config, plotData,dictPlot, bosonName, targetvariable, resultData, dictResult, minrange=42,maxrange=0, stepwidth=0, ptmin=0, ptmax=0):
  
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
    print('Response Perp %s versus %s'%(bosonName,targetvariable))
    
    #YValues 
    ignoredEntries = 0
    for index in range(0,XRange.shape[0]-1):

	AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]

	currentDistri = AlternativeDistri[:,dictPlot[bosonName]]
	"""
	if AlternativeDistri.shape[0] == 0:
	    YMean[index] = 0
	    YStd[index] = 0
	else:
	"""
	YMean[index] = currentDistri.mean()
	YStd[index] = currentDistri.std()
	if index < 12:
	    plt.clf()
            num_bins = 50
	    n, bins, patches = plt.hist(currentDistri, num_bins, normed=1, facecolor='green', alpha=0.5)
	    y = mlab.normpdf(bins, currentDistri.mean(), currentDistri.std())
	    plt.xlabel('%s at %f, mean: %f'%(targetvariable,(XRange[index+1]+XRange[index])/2,currentDistri.mean()))
	    plt.ylabel('MET Boson PT_Perp')
	    plt.plot(bins, y, 'r--')
	    if ptmax == 0:
		plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/ResponsePerp_%s_vs_%s_%i.png' %(bosonName,targetvariable, index)))
	    else:
		plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/ResponsePerp(%i<Pt<%i)_%s_vs_%s_%i.png' %(ptmin,ptmax,bosonName,targetvariable, index)))
    plt.clf()
    plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o')

    plt.xlabel(targetvariable)
    if ptmax == 0:
	plt.savefig(config['outputDir'] + "ControlPlots/ResponsePerp_%s_vs_%s.png" %(bosonName,targetvariable))
    else:
	plt.savefig(config['outputDir'] + "ControlPlots/ResponsePerp(%i<Pt<%i)_%s_vs_%s.png" %(ptmin,ptmax,bosonName,targetvariable))
    plt.clf()
    plt.figure(7)
    plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o',label=bosonName)
    plt.figure(0)
    

    if bosonName == 'LongZCorrectedRecoil_PerpZ':
    #if bosonName == 'recoilslimmedMETs_LongZ':
	if ptmax > config[config['activePlotting']]['BosonCut']:
	    Upper = 'Max'
	else:
	    Upper = 'Cut'
	    
	if ptmin < config[config['activePlotting']]['BosonCut']:
	    Lower = 'Min'
	else:
	    Lower = 'Cut'
      
      
	if targetvariable == 'Boson_Pt':
	    resultvalue = ((YMean[np.isfinite(YMean)]-1)**2).sum()/YMean[np.isfinite(YMean)].shape[0]
	    dictResult["ResponsePerp_Chi_%s_MintoMax"%(targetvariable)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
	    XRangeTemp = XRange[:-1]
	    
	    
	    YMeanCut = YMean[config[config['activePlotting']]['BosonCut']<XRangeTemp[:]]
	    resultvalue = ((YMeanCut[np.isfinite(YMeanCut)]-1)**2).sum()/YMeanCut[np.isfinite(YMeanCut)].shape[0]
	    dictResult["ResponsePerp_Chi_%s_CuttoMax"%(targetvariable)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
	elif targetvariable == 'NVertex':
	    
	    resultvalue = ((YMean[np.isfinite(YMean)]-1)**2).sum()/YMean[np.isfinite(YMean)].shape[0]
	    dictResult["ResponsePerp_Chi_%s_%sto%s"%(targetvariable,Lower,Upper)] = resultData.shape[1]
	    resultData = np.hstack((resultData, np.array(resultvalue.reshape(1,1))))
    
    return resultData, dictResult

    

def make_ControlPlots(config, plotData,dictPlot, bosonName, targetvariable, resultData, dictResult, minrange=42,maxrange=0, stepwidth=0, ptmin=0,ptmax=0):

    bosonNameLong = bosonName + '_LongZ'
    bosonNamePerp = bosonName + '_PerpZ'
    maxrange += stepwidth
    if not os.path.exists((config['outputDir'] + 'ControlPlots/SingleDistributions/')):
	os.makedirs((config['outputDir'] + 'ControlPlots/SingleDistributions/'))
    resultData, dictResult, XRange, YVariance = make_ResolutionPlot(config, plotData, dictPlot, bosonNameLong, targetvariable, resultData, dictResult, minrange,maxrange,stepwidth, ptmin, ptmax) 
    resultData, dictResult, YResponse = make_ResponsePlot(config, plotData, dictPlot, bosonNameLong, targetvariable, resultData, dictResult, minrange,maxrange,stepwidth, ptmin, ptmax)
    resultData, dictResult = make_ResponseCorrectedPlot(config, XRange, YVariance, YResponse, bosonNameLong, targetvariable, resultData, dictResult, minrange,maxrange, stepwidth, ptmin, ptmax)
    resultData, dictResult = make_ResolutionPerpPlot(config, plotData, dictPlot, bosonNamePerp, targetvariable, resultData, dictResult, minrange,maxrange,stepwidth, ptmin, ptmax) 
    resultData, dictResult = make_ResponsePerpPlot(config, plotData, dictPlot, bosonNamePerp, targetvariable, resultData, dictResult, minrange,maxrange,stepwidth, ptmin, ptmax)
    
    if bosonName == "LongZCorrectedRecoil" and "recoMetOnGenMetProjectionPar" in dictPlot and not plotData[dictPlot["recoMetOnGenMetProjectionPar"]].mean() == -999:
        make_METResponsePlot(config, plotData, dictPlot, "recoMetOnGenMetProjectionPar", targetvariable, resultData, dictResult, minrange,maxrange,stepwidth, ptmin, ptmax)
        make_METResolutionPlot(config, plotData, dictPlot, "recoMetOnGenMetProjectionPar", targetvariable, resultData, dictResult, minrange,maxrange,stepwidth, ptmin, ptmax) 
    
    if bosonName == "recoilslimmedMETs" and "recoPfMetOnGenMetProjectionPar" in dictPlot and not plotData[dictPlot["recoPfMetOnGenMetProjectionPar"]].mean() == -999:
        make_METResponsePlot(config, plotData, dictPlot, "recoPfMetOnGenMetProjectionPar", targetvariable, resultData, dictResult, minrange,maxrange,stepwidth, ptmin, ptmax)
        make_METResolutionPlot(config, plotData, dictPlot, "recoPfMetOnGenMetProjectionPar", targetvariable, resultData, dictResult, minrange,maxrange,stepwidth, ptmin, ptmax) 
    
    
    
    return resultData, dictResult

    
def make_MoreBDTPlots(config, plotData, dictPlot):
    
    plt.clf()
    plt.hist2d(plotData[:,dictPlot['Boson_Pt']], plotData[:,dictPlot['flatPtWeight']],bins = 80, norm=LogNorm())
    plt.xlabel('Boson Pt')
    plt.ylabel('flat Pt weight')
    plt.savefig(config['outputDir'] + "/ControlPlots/BosonOverWeights.png")
    
    num_bins = 50
    
    plt.clf()
    if 'fileName' in dictPlot:
	for i in range (0,9):
	    if i in plotData[:,dictPlot['fileName']]:
		currentDistri = plotData[i==plotData[:,dictPlot['fileName']],dictPlot['fileName']]
		n, bins, patches = plt.hist(currentDistri, num_bins, range=[0,plotData[:,dictPlot['Boson_Pt']].max()], alpha=0.5, label=('File %i'%i))
	plt.xlabel('Boson Pt')
	plt.ylabel('Entries')
	plt.savefig(config['outputDir'] + "/ControlPlots/BosonPtSpectrum.png")
	plt.clf()
      
    histBosonPtLow = plotData[plotData[:,dictPlot['Boson_Pt']]<30,dictPlot['Boson_Pt']]
    
    n, bins, patches = plt.hist(histBosonPtLow, num_bins, range=[0,30],alpha=0.5)
    plt.xlabel('Boson Pt')
    plt.ylabel('Entries')
    plt.savefig(config['outputDir'] + "/ControlPlots/LowBosonPtSpectrum.png")
    plt.clf()

    borders = [0,10,20,30,40,50,60,70,80,90,100]
    for index,min in enumerate(borders):
	histname = 'targetRecoilFromBDT'
	slicedData = plotData[min<=plotData[:,dictPlot['Boson_Pt']],:]
	slicedData = slicedData[min+10>=slicedData[:,dictPlot['Boson_Pt']],:]
	histDataTargetBDT = slicedData[:,dictPlot[histname]]

	histDataOutputBDT = slicedData[:,dictPlot['LongZCorrectedRecoil_LongZ']]/slicedData[:,dictPlot['PhiCorrectedRecoil_LongZ']]
		
	histDataVarianceBDT = histDataOutputBDT - histDataTargetBDT
	
	histDataAtOnce = [histDataTargetBDT,histDataOutputBDT]
		
	names = ['Target scale factor, mean: %f'%histDataTargetBDT.mean(),'BDT predicted scale factor, mean %f'%histDataOutputBDT.mean()]
	n, bins, patches = plt.hist(histDataAtOnce, num_bins, range=[-2, 4], alpha=0.5, label=names)
	plt.legend(loc='upper right')
	plt.xlabel('Comparison scalefactor BDT Output and Target. %i GeV < Boson Pt < %i'%(min,min+10))
	plt.ylabel('Entries')
	plt.savefig(config['outputDir'] + "/ControlPlots/BDT_Scalefactor_OutputAndTarget%ito%i.png"%(min,min+10))
	plt.clf()

    
    histVarJet0AndBosonPtCut = np.abs(plotData[30<plotData[:,dictPlot['Boson_Pt']],dictPlot['Boson_Phi']] - plotData[30<plotData[:,dictPlot['Boson_Pt']],dictPlot['Jet0_Phi']])

    n, bins, patches = plt.hist(histVarJet0AndBosonPtCut, num_bins,alpha=0.5)
    plt.xlabel('Var Phi (Jet0,BosonPt), BosonPt > 30')
    plt.ylabel('Entries')
    plt.savefig(config['outputDir'] + "/ControlPlots/VarJet0AndBosonPtCut30.png")
    plt.clf()
    
    histVarJet1AndBosonPtCut = np.abs(plotData[30<plotData[:,dictPlot['Boson_Pt']],dictPlot['Boson_Phi']] - plotData[30<plotData[:,dictPlot['Boson_Pt']],dictPlot['Jet1_Phi']])
    

    n, bins, patches = plt.hist(histVarJet1AndBosonPtCut, num_bins,alpha=0.5)
    plt.xlabel('Var Phi (Jet1,BosonPt), BosonPt > 30')
    plt.ylabel('Entries')
    plt.savefig(config['outputDir'] + "/ControlPlots/VarJet1AndBosonPtCut30.png")
    plt.clf()
    
    histVarJet0AndBosonCleaned = np.abs(plotData[10<plotData[:,dictPlot['Jet0_Pt']],dictPlot['Boson_Phi']] - plotData[10<plotData[:,dictPlot['Jet0_Pt']],dictPlot['Jet0_Phi']])

    n, bins, patches = plt.hist(histVarJet0AndBosonCleaned, num_bins,alpha=0.5)
    plt.xlabel('Var Phi (Jet0,BosonPt), JetPt > 10')
    plt.ylabel('Entries')
    plt.savefig(config['outputDir'] + "/ControlPlots/VarJet0AndBosonCleaned10.png")
    plt.clf()
    
    histVarJet1AndBosonCleaned = np.abs(plotData[10<plotData[:,dictPlot['Jet1_Pt']],dictPlot['Boson_Phi']] - plotData[10<plotData[:,dictPlot['Jet1_Pt']],dictPlot['Jet1_Phi']])
    

    n, bins, patches = plt.hist(histVarJet1AndBosonCleaned, num_bins,alpha=0.5)
    plt.xlabel('Var Phi (Jet1,BosonPt), JetPt > 10')
    plt.ylabel('Entries')
    plt.savefig(config['outputDir'] + "/ControlPlots/VarJet1AndBosonCleaned10.png")
    plt.clf()
    
    
    histVarJet0AndBoson = np.abs(plotData[:,dictPlot['Boson_Phi']] - plotData[:,dictPlot['Jet0_Phi']])

    n, bins, patches = plt.hist(histVarJet0AndBoson, num_bins,alpha=0.5)
    plt.xlabel('Var Phi (Jet0,BosonPt)')
    plt.ylabel('Entries')
    plt.savefig(config['outputDir'] + "/ControlPlots/VarJet0AndBoson.png")
    plt.clf()
    
    histVarJet1AndBoson = np.abs(plotData[:,dictPlot['Boson_Phi']] - plotData[:,dictPlot['Jet1_Phi']])
    

    n, bins, patches = plt.hist(histVarJet1AndBoson, num_bins,alpha=0.5)
    plt.xlabel('Var Phi (Jet1,BosonPt)')
    plt.ylabel('Entries')
    plt.savefig(config['outputDir'] + "/ControlPlots/VarJet1AndBoson.png")
    plt.clf()
    
    histJet0PtSmall = plotData[:,dictPlot['Jet0_Pt']]
    
    n, bins, patches = plt.hist(histJet0PtSmall, num_bins, range=[0, 100],alpha=0.5)
    plt.xlabel('Jet 0 Pt')
    plt.ylabel('Entries')
    plt.savefig(config['outputDir'] + "/PlotVariables/Jet0PtSmall.png")
    plt.clf()
    
    
    
    if 'select' in dictPlot:
	muData = plotData[plotData[:,dictPlot['select']]==1]
	eData = plotData[plotData[:,dictPlot['select']]==2]
	
	plotVariable = 'Boson_Pt'
	
	histData = [muData[:,dictPlot[plotVariable]],eData[:,dictPlot[plotVariable]]]
	names = ['%s - Z to mumu'%plotVariable, '%s - Z to ee'%plotVariable]
	n, bins, patches = plt.hist(histData, num_bins, range=[0, 300], normed=1, alpha=0.5, label=names)
	plt.legend(loc='upper right')
	plt.xlabel('Boson Pt in GeV')
	plt.ylabel('Entries')
	plt.savefig(config['outputDir'] + "/ControlPlots/BDT_MuEComparison_%s"%plotVariable)
	plt.clf()

	plotVariable = 'NVertex'
	num_bins = 35
	histData = [muData[:,dictPlot[plotVariable]],eData[:,dictPlot[plotVariable]]]
	names = ['%s - Z to mumu'%plotVariable, '%s - Z to ee'%plotVariable]
	n, bins, patches = plt.hist(histData, num_bins, range=[5, 40], normed=1, alpha=0.5, label=names)
	plt.legend(loc='upper right')
	plt.xlabel('Boson Pt in GeV')
	plt.ylabel('Entries')
	plt.savefig(config['outputDir'] + "/ControlPlots/BDT_MuEComparison_%s"%plotVariable)
	plt.clf()
	
	
	plotVariable = 'NCleanedJets'
	
	num_bins = 15
	histData = [muData[:,dictPlot[plotVariable]],eData[:,dictPlot[plotVariable]]]
	names = ['%s - Z to mumu'%plotVariable, '%s - Z to ee'%plotVariable]
	n, bins, patches = plt.hist(histData, num_bins, range=[0, 15], normed=1, alpha=0.5, label=names)
	plt.legend(loc='upper right')
	plt.xlabel('Boson Pt in GeV')
	plt.ylabel('Entries')
	plt.savefig(config['outputDir'] + "/ControlPlots/BDT_MuEComparison_%s"%plotVariable)
	plt.clf()
    
    


def make_PtSpectrumPlot(config, plotData, dictPlot, maxBosonPt=0, stepwidth=0):

  
    if maxBosonPt==0:
	maxBosonPt=plotData[:,dictPlot['Boson_Pt']].max()
    if stepwidth==0:
	stepwidth = maxBosonPt/100
    
    XRange = np.arange(0,maxBosonPt,stepwidth)
    YSum = np.zeros((XRange.shape[0]-1,1))
    for index in range(0,XRange.shape[0]-1):
	AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot['Boson_Pt']]) & (XRange[index+1]>plotData[:,dictPlot['Boson_Pt']]) & (plotData[:,dictPlot['nCombinations']]==1)]
	sumEntries = AlternativeDistri[:,dictPlot['flatPtWeight']].sum()
	YSum[index] = sumEntries
 
    plt.clf()
    plt.plot(XRange[:-1],YSum[:],'o')
    plt.xlabel('Boson Pt')
    plt.ylabel('Weighted Boson Pt')
    plt.savefig(config['outputDir'] + "WeightedBosonPt.png")
    
    plt.clf()
    
    weightPt = plotData[:,dictPlot['flatPtWeight']]
    
    num_bins = 50
    n, bins, patches = plt.hist(plotData[:,dictPlot['Boson_Pt']], num_bins, facecolor='green', alpha=0.5, weights=weightPt)
    plt.savefig(config['outputDir'] + "WeightedBosonPtHist.png")
    plt.clf()
    
    if 'flatPtWeightBDT' in dictPlot:
	n, bins, patches = plt.hist(plotData[:,dictPlot['Boson_Pt']], num_bins, facecolor='green', alpha=0.5, weights=plotData[:,dictPlot['flatPtWeightBDT']])
	plt.savefig(config['outputDir'] + "WeightedBosonPtHistBDT.png")
	plt.clf()
    
 
 

def make_PhiVariancePlot(config, plotData, dictPlot, targetvariable, ptmin, ptmax, xlabelname = ''):
  
    if xlabelname == '':
	xlabelname = targetvariable
    num_bins = 50
    histDataPhi = plotData[:,dictPlot['Boson_Phi']] + math.pi - plotData[:,dictPlot[targetvariable]]
    print(targetvariable,' phi shape: ',histDataPhi.shape)
    for event in range(0,histDataPhi.shape[0]):
	if histDataPhi[event] > math.pi:
	    histDataPhi[event] -= 2*math.pi
	if histDataPhi[event] < -math.pi:
	    histDataPhi[event] += 2*math.pi
    MSE = (histDataPhi**2).mean()
    n, bins, patches = plt.hist(histDataPhi, num_bins, facecolor='green', alpha=0.5)
    plt.xlabel('Variance %s from true Boson Phi (%i < Boson Pt < %i)GeV. MSE: %f'%(xlabelname, ptmin, ptmax, MSE))
    plt.ylabel('Entries')
    plt.savefig(config['outputDir'] + "PhiVariance%s_%ito%iPt.png"%(xlabelname,ptmin,ptmax))
    plt.clf()
    print('MSE %s (%i<Pt<%i): '%(xlabelname,ptmin,ptmax),MSE)    

    # normal distribution center at x=0 and y=5
    plt.hist2d(plotData[:,dictPlot['Boson_Phi']], histDataPhi,bins = 80, norm=LogNorm())
    #plt.ylim([-0.25,0.25])
    plt.xlabel('Boson Phi (%i < Boson Pt < %i)GeV'%(ptmin, ptmax))
    plt.ylabel('Variance of (Prediction-Target) %s'%xlabelname)
    plt.savefig(config['outputDir'] + "Variance2D_%s(%i<Pt<%i).png"%(xlabelname,ptmin,ptmax))
    plt.clf()

    
    
    
def plot_results(config, plotData, dictPlot):
    
    #plotData = plotData[0==plotData[:,dictPlot['NCleanedJets']],:]
    
    resultData = np.empty(shape=[1,0]).astype(np.float32)
    dictResult = {}
    
    
    
    
    
    
    
    plotconfig = config[config['activePlotting']]
    
    num_bins = 50
  
    if not os.path.exists(config['outputDir']):
	os.makedirs(config['outputDir'])
	
    
    #if not os.path.exists('../plots/PlotVariables'):
	#os.makedirs('../plots/PlotVariables')
    
    if plotconfig['plotPlotVariables']:
	if plotconfig['splitNCleanedJets']:
            if 'NCleanedJets' in dictPlot:
                plotVariableData = plotData[0<plotData[:,dictPlot['NCleanedJets']],:]
                outputdir = config['outputDir'] + 'PlotVariables/1CleanedJet/'
                for variable in dictPlot:
                    make_Plot(variable,plotVariableData,dictPlot,outputdir)
                    
                plotVariableData = plotData[1<plotData[:,dictPlot['NCleanedJets']],:]
                outputdir = config['outputDir'] + 'PlotVariables/2CleanedJet/'
                for variable in dictPlot:
                    make_Plot(variable,plotVariableData,dictPlot,outputdir)
                    
                plotVariableData = plotData[0==plotData[:,dictPlot['NCleanedJets']],:]
                outputdir = config['outputDir'] + 'PlotVariables/UncleanedJet/'
                for variable in dictPlot:
                    make_Plot(variable,plotVariableData,dictPlot,outputdir)
		
	outputdir = config['outputDir'] + 'PlotVariables/'
	for variable in dictPlot:
	    make_Plot(variable,plotData,dictPlot,outputdir)
      
    #Boson Pt
    #comparisonMinus.xlabel('Boson_Pt')
    #comparisonOver.xlabel('Boson_Pt')
    #comparisonMinus.ylabel('|Boson_Pt-Prediction|')
    #comparisonOver.ylabel('Prediction/Boson_Pt')
    
    bosonmin = [0,0,plotconfig['BosonCut']]
    bosonmax = [plotData[:,dictPlot['Boson_Pt']].max(),plotconfig['BosonCut'],plotData[:,dictPlot['Boson_Pt']].max()]
    
    #BDT Performance Plots 
#Phi ControlPlots
    for i, min in enumerate(bosonmin):
        slicedData = plotData[min<=plotData[:,dictPlot['Boson_Pt']],:]
        slicedData = slicedData[bosonmax[i]>=slicedData[:,dictPlot['Boson_Pt']],:]
        if plotconfig['plotBDTPerformance']:
            if plotconfig['plotAdditionalBDTPlots']:
                histname = 'targetRecoilFromBDT'

                histDataTargetBDT = slicedData[:,dictPlot[histname]]

                histDataOutputBDT = slicedData[:,dictPlot['LongZCorrectedRecoil_LongZ']]/slicedData[:,dictPlot['PhiCorrectedRecoil_LongZ']]
                
                histDataVarianceBDT = histDataOutputBDT - histDataTargetBDT
                
                n, bins, patches = plt.hist(histDataVarianceBDT, num_bins, range=[-2, 4], facecolor='green', alpha=0.5)
                plt.xlabel('Scalefactor variance (BDT-Target). (%i  < Boson Pt < %i)GeV, mean: %f'%(min,bosonmax[i],histDataVarianceBDT.mean()))
                plt.ylabel('Entries')
                plt.savefig(config['outputDir'] + "BDT_Scalefactor_VarianceOutputandTarget%ito%i.png"%(min,bosonmax[i]))
                plt.clf()
                
                histDataAtOnce = [histDataTargetBDT,histDataOutputBDT]
                
                names = ['Target scale factor, mean: %f'%histDataTargetBDT.mean(),'BDT predicted scale factor, mean %f'%histDataOutputBDT.mean()]
                n, bins, patches = plt.hist(histDataAtOnce, num_bins, range=[-2, 4], alpha=0.5, label=names)
                plt.legend(loc='upper right')
                plt.xlabel('Comparison scalefactor BDT Output and Target. %i GeV < Boson Pt < %i'%(min,bosonmax[i]))
                plt.ylabel('Entries')
                plt.savefig(config['outputDir'] + "BDT_Scalefactor_OutputAndTarget%ito%i.png"%(min,bosonmax[i]))
                plt.clf()
    
                make_PhiVariancePlot(config, slicedData,dictPlot,'PhiCorrectedRecoil_Phi', min, bosonmax[i], 'BDT Phi')
    
    
            #BDT
            if 'LongZCorrectedRecoil_LongZ' in dictPlot:
                resultData, dictResult = make_ControlPlots(config, slicedData, dictPlot, 'LongZCorrectedRecoil', 'NVertex',resultData, dictResult, 5,40,5,min,bosonmax[i])
            
            if 'PhiCorrectedRecoil_LongZ' in dictPlot:
                resultData, dictResult = make_ControlPlots(config, slicedData, dictPlot, 'PhiCorrectedRecoil', 'NVertex',resultData, dictResult, 5,40,5,min,bosonmax[i])
            
        
        
        
        #slimmedMet
        if 'recoilslimmedMETs_LongZ' in dictPlot:
            resultData, dictResult = make_ControlPlots(config, slicedData, dictPlot, 'recoilslimmedMETs', 'NVertex',resultData, dictResult, 5,40,5,min,bosonmax[i])
        
        #Puppi-Met
        if 'recoilslimmedMETsPuppi_LongZ' in dictPlot:
            resultData, dictResult = make_ControlPlots(config, slicedData, dictPlot, 'recoilslimmedMETsPuppi', 'NVertex',resultData, dictResult, 5,40,5,min,bosonmax[i])
        
        plt.figure(4)
        legend = plt.legend(loc='lower right', shadow=True)
        #legend.get_frame().set_alpha(0.5)
        plt.xlabel('N Vertex (%i < Boson Pt < %i)'%(min,bosonmax[i]))
        plt.ylabel('(MET Boson PT_Long)/(True Boson Pt)')
        plt.savefig(config['outputDir'] + 'Response_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()
        plt.figure(5)
        plt.xlabel('N Vertex (%i < Boson Pt < %i)'%(min,bosonmax[i]))
        plt.ylabel('(MET Boson PT_Long) - (True Boson Pt)')
        legend = plt.legend(loc='lower right', shadow=True)
        #legend.get_frame().set_alpha(0.5)
        plt.savefig(config['outputDir'] + 'Resolution_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()
        plt.figure(6)
        plt.xlabel('N Vertex (%i < Boson Pt < %i)'%(min,bosonmax[i]))
        plt.ylabel('Resolution / Response')
        legend = plt.legend(loc='lower right', shadow=True)
        #legend.get_frame().set_alpha(0.5)
        plt.savefig(config['outputDir'] + 'ResponseCorrected_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()
        plt.figure(7)
        legend = plt.legend(loc='lower right', shadow=True)
        #legend.get_frame().set_alpha(0.5)
        plt.xlabel('N Vertex (%i < Boson Pt < %i)'%(min,bosonmax[i]))
        plt.ylabel('MET Boson PT_Perp')
        plt.savefig(config['outputDir'] + 'ResponsePerp_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()
        plt.figure(8)
        plt.xlabel('N Vertex (%i < Boson Pt < %i)'%(min,bosonmax[i]))
        plt.ylabel('MET Boson PT_Perp')
        legend = plt.legend(loc='lower right', shadow=True)
        #legend.get_frame().set_alpha(0.5)
        plt.savefig(config['outputDir'] + 'ResolutionPerp_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()
        plt.figure(9)
        plt.xlabel('N Vertex (%i < Boson Pt < %i)'%(min,bosonmax[i]))
        plt.ylabel('MET/genMET')
        legend = plt.legend(loc='lower right', shadow=True)
        #legend.get_frame().set_alpha(0.5)
        plt.savefig(config['outputDir'] + 'METResponse_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()
        plt.figure(10)
        plt.xlabel('N Vertex (%i < Boson Pt < %i)'%(min,bosonmax[i]))
        plt.ylabel('MET - genMET')
        legend = plt.legend(loc='lower right', shadow=True)
        #legend.get_frame().set_alpha(0.5)
        plt.savefig(config['outputDir'] + 'Resolution_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()
        plt.figure(0)
        

        if plotconfig['plotAdditionalBDTPlots']:
            make_PhiVariancePlot(config, slicedData,dictPlot,'recoilslimmedMETs_Phi', min, bosonmax[i], 'PF Phi')
        
            make_PhiVariancePlot(config, slicedData,dictPlot,'recoilslimmedMETsPuppi_Phi', min, bosonmax[i],'PUPPI Phi')
	  
    
    
    #Boson PT
    

    if 'recoilslimmedMETs_LongZ' in dictPlot:
        resultData, dictResult = make_ControlPlots(config, plotData, dictPlot, 'recoilslimmedMETs',                       'Boson_Pt',resultData, dictResult, 10,200,10,0,0)
    
    if plotconfig['plotBDTPerformance']:
        if 'LongZCorrectedRecoil_LongZ' in dictPlot:
            resultData, dictResult = make_ControlPlots(config, plotData, dictPlot, 'LongZCorrectedRecoil', 'Boson_Pt',resultData, dictResult, 10,200,10,0,0)
        if 'PhiCorrectedRecoil_LongZ' in dictPlot:
            resultData, dictResult = make_ControlPlots(config, plotData, dictPlot, 'PhiCorrectedRecoil', 'Boson_Pt',resultData, dictResult, 10,200,10,0,0)
        if plotconfig['plotAdditionalBDTPlots']:
            make_MoreBDTPlots(config, plotData, dictPlot)

    
    
    if 'recoilslimmedMETsPuppi_LongZ' in dictPlot:
        resultData, dictResult = make_ControlPlots(config, plotData, dictPlot, 'recoilslimmedMETsPuppi', 'Boson_Pt',resultData, dictResult, 10,200,10,0,0)


    plt.figure(4)
    legend = plt.legend(loc='lower right', shadow=True)
    #legend.get_frame().set_alpha(0.5)
    plt.xlabel('Boson Pt')
    plt.ylabel('(MET Boson PT_Long)/(True Boson Pt)')
    plt.savefig(config['outputDir'] + 'Response_vs_BosonPt')
    plt.clf()
    plt.figure(5)
    plt.xlabel('Boson Pt')
    plt.ylabel('(MET Boson PT_Long) - (True Boson Pt)')
    legend = plt.legend(loc='lower right', shadow=True)
    #legend.get_frame().set_alpha(0.5)
    plt.savefig(config['outputDir'] + 'Resolution_vs_BosonPt')
    plt.clf()
    plt.figure(6)
    plt.xlabel('Boson Pt')
    plt.ylabel('Resolution / Response')
    legend = plt.legend(loc='lower right', shadow=True)
    #legend.get_frame().set_alpha(0.5)
    plt.savefig(config['outputDir'] + 'ResponseCorrected_vs_BosonPt')
    plt.clf()
    plt.figure(7)
    legend = plt.legend(loc='lower right', shadow=True)
    #legend.get_frame().set_alpha(0.5)
    plt.xlabel('Boson Pt')
    plt.ylabel('MET Boson PT_Perp')
    plt.savefig(config['outputDir'] + 'ResponsePerp_vs_BosonPt')
    plt.clf()
    plt.figure(8)
    plt.xlabel('Boson Pt')
    plt.ylabel('MET Boson PT_Perp')
    legend = plt.legend(loc='lower right', shadow=True)
    #legend.get_frame().set_alpha(0.5)
    plt.savefig(config['outputDir'] + 'ResolutionPerp_vs_BosonPt')
    plt.clf()
    plt.figure(9)
    legend = plt.legend(loc='lower right', shadow=True)
    #legend.get_frame().set_alpha(0.5)
    plt.xlabel('Boson Pt')
    plt.ylabel('MET/genMET')
    plt.savefig(config['outputDir'] + 'METResponse_vs_BosonPt')
    plt.clf()
    plt.figure(10)
    plt.xlabel('Boson Pt')
    plt.ylabel('MET - genMET')
    legend = plt.legend(loc='lower right', shadow=True)
    #legend.get_frame().set_alpha(0.5)
    plt.savefig(config['outputDir'] + 'METResolution_vs_BosonPt')
    plt.clf()
    plt.figure(0)
	    
    
    print('resultDatashape: ',resultData.shape)
    for name in dictResult:
	print(name)
	print(resultData[0,dictResult[name]])
    
    
    
    #make Boson Pt Spectrum with weights
    if plotconfig['plotAdditionalBDTPlots']:
        make_PtSpectrumPlot(config, plotData, dictPlot, 650, 5)

    print('Plots created')
    
    return True

    



def main(config):
    

    # Load the dataset
    print("Loading data...")
    plotData, dictPlot = load_datasetcsv(config)
    

    print("plotDatashape", plotData.shape)
    plot_results(config, plotData, dictPlot)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make MVAMet control plots.')
    parser.add_argument('-config', '--inputconfig', default='../configs/config.json', help='Path to configurations file')
    parser.add_argument('-i', '--inputfile', default='NaN', help='[optional] Inputfile from which to create the plots from')
    args = parser.parse_args()
    print('Used Configfile: ',args.inputconfig)



    with open(args.inputconfig, 'r') as f:
	config = json.load(f)
    
    if not args.inputfile == 'NaN':
        config['inputFile'] = args.inputfile
        
    print('Inputfile: %s'%config['inputFile'])
    
    main(config)


