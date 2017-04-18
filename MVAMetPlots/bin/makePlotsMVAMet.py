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
from scipy.interpolate import splrep, sproot, splev
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
def get_Quantity(data, dictPlot, quantity, method = 'Trunc'):
    if method == 'Empiric':
	if data.shape[0] == 0:
	    return 0
	else:
	    if quantity == 'Mean':
		return data.mean()
	    elif quantity == 'Std':
		return data.std()
	    else:
		return 0
    else:
	centralData = data[((data.mean()-4*data.std())<data[:]) & (data[:] <(data.mean()+4*data.std()))]
	if centralData.shape[0] == 0:
	    return 0
	else:
	    if method == 'Trunc':
		if quantity == 'Mean':
		    return centralData.mean()
		elif quantity == 'Std':
		    return centralData.std()
		else:
		    return 0
	    else:
		num_bins = 50
		n, bins, patches = plt.hist(centralData, num_bins, facecolor='green', alpha=0.5)
		XCenter = (bins[:-1] + bins[1:])/2
		if method == 'Fit':
		    p0 = [1., 0., 1.]
		    try:
			coeff, var_matrix = curve_fit(gauss,XCenter,n,p0=p0)
		    except:
			coeff =[0,0,0]
		    if quantity == 'Mean':
			return coeff[1]
		    elif quantity == 'Std':
			return abs(coeff[2])
		    else:
			return 0
		elif method == 'FWHM':
		    FWHM_mean, FWHM_std = fwhm(XCenter,n)
		    if quantity == 'Mean':
			return FWHM_mean
		    elif quantity == 'Std':
			return FWHM_std/2.355
		    else:
			return 0
		elif method == 'FWHMSpline':
                    FWHM_mean, FWHM_std = fwhm(XCenter,n,'Spline')
                    if quantity == 'Mean':
                        return FWHM_mean
                    elif quantity == 'Std':
                        return FWHM_std/2.355
                    else:
                        return 0
		else:
		    return 0

def getNextAlphaPoint(dictPlot, alphaData, currentPoint):
	AlphaVar = alphaData[:,dictPlot['Jet1_Pt']]/alphaData[:,dictPlot['Boson_Pt']]
	AlphaVar = AlphaVar[1>AlphaVar[:]]
	if currentPoint == 0:
		return np.min(AlphaVar)
	else:
		entriesAlpha = alphaData.shape[0]
		deltaEntries = entriesAlpha/10
		i = 1
		while alphaData[(alphaData[:,dictPlot['Jet1_Pt']]/alphaData[:,dictPlot['Boson_Pt']]>currentPoint) &(alphaData[:,dictPlot['Jet1_Pt']]/alphaData[:,dictPlot['Boson_Pt']]<(currentPoint+0.01*i))].shape[0] < deltaEntries and (currentPoint + 0.01*i)<= 1.:
			i += 1
		return currentPoint + i * 0.01


def alphaFit(dictPlot, alphaDataIn, bosonName='recoilslimmedMETs_Pt', mode='Response', saveName = 'dummy.png', method = 'AlphaInclInt'):
	XRange = np.zeros(9)

    	YMean = np.zeros(9)
    	YStd = np.zeros(9)

	alphaData = alphaDataIn[(0<alphaDataIn[:,dictPlot['Jet1_Pt']])]
	try:
	    minAlpha = getNextAlphaPoint(dictPlot, alphaData, 0)
	    maxAlpha = getNextAlphaPoint(dictPlot, alphaData, minAlpha)
	    for index in range(9):
		    if method == 'AlphaInclInt':
			alphaDataLoop = alphaDataIn[(alphaDataIn[:,dictPlot['Jet1_Pt']]/alphaDataIn[:,dictPlot['Boson_Pt']]<maxAlpha)]
		    elif method == 'AlphaExclInt':
			alphaDataLoop = alphaData[(alphaData[:,dictPlot['Jet1_Pt']]/alphaData[:,dictPlot['Boson_Pt']]<maxAlpha)]
		    elif method == 'AlphaExcl':
			alphaDataLoop = alphaData[(alphaData[:,dictPlot['Jet1_Pt']]/alphaData[:,dictPlot['Boson_Pt']]>minAlpha) &(alphaData[:,dictPlot['Jet1_Pt']]/alphaData[:,dictPlot['Boson_Pt']]<maxAlpha)]

		    if mode == 'Response':
			    currentDistri = -alphaDataLoop[:,dictPlot[bosonName]]/alphaDataLoop[:,dictPlot['Boson_Pt']]
		    elif mode == 'Resolution':
			    currentDistri = alphaDataLoop[:,dictPlot[bosonName]]+alphaDataLoop[:,dictPlot['Boson_Pt']]
		    #fitDistri = currentDistri[((currentDistri.mean()-4*currentDistri.std())<currentDistri[:]) & (currentDistri[:] <(currentDistri.mean()+4*currentDistri.std()))]
		    #XRange[index] = (maxAlpha+minAlpha)/2
		    XRange[index] = maxAlpha
		    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', method = 'Trunc')
		    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', method = 'Trunc')
		    minAlpha = maxAlpha
		    maxAlpha = getNextAlphaPoint(dictPlot, alphaData, minAlpha)
	except:
	    print('Stefan nervt')

	p0 = [0.,1.]

	plt.clf()
	try:
		coeffMean, var_matrix = curve_fit(linear,XRange.transpose(),YMean.transpose(),p0=p0)
	except:
		coeffMean = [0.,0.]
	try:
		coeffStd, var_matrix = curve_fit(linear,XRange.transpose(),YStd.transpose(),p0=p0)
	except:
		coeffStd = [0.,0.]

	if mode == 'Response':
		plt.plot(XRange,YMean,'o')
		y = linear(XRange,*coeffMean)
		plt.plot(XRange,y)
		plt.ylabel(r'$<U_\| / p_t^Z>$'' in GeV',fontsize = 17)
		plt.text(0.60*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.30*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$\mathrm{Events\,per\,Fit} = %i$''\n'r'$y = a \cdot x + b$''\n'r'$\mathrm{a} = %.2f$''\n'r'$\mathrm{b} = %.2f$''\n'%(alphaData.shape[0]/10,coeffMean[0], coeffMean[1]),color = 'k',fontsize=16)
	elif mode == 'Resolution':
		plt.plot(XRange,YStd,'o')
		y = linear(XRange,*coeffStd)
		plt.plot(XRange,y)
		plt.ylabel(r'$\sigma(<U_\| - p_t^Z>)$'' in GeV',fontsize = 17)
		plt.text(0.60*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.30*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$\mathrm{Events\,per\,Fit} = %i$''\n'r'$y = a \cdot x + b$''\n'r'$\mathrm{a} = %.2f$''\n'r'$\mathrm{b} = %.2f$''\n'%(alphaData.shape[0]/10,coeffStd[0], coeffStd[1]),color = 'k',fontsize=16)
	plt.xlabel(r'$\alpha = p_{t}^{Jet1}/p_t^Z}$',fontsize = 17)
	if not saveName == 'dummy.png':
		plt.savefig(saveName)
	plt.clf()

	return coeffMean[1], coeffStd[1]


def fwhm(x, y, method = 'Linear',k=10):
    """
    Determine full-with-half-maximum of a peaked set of points, x and y.

    Assumes that there is only one peak present in the datasset.  The function
    uses a spline interpolation of order k.
    """



    half_max = np.max(y)/2.0

    if method == 'Linear':
	#version B / linear interpolation
	roots = []
	for index in range(len(x)-1):
	    if y[index] <= half_max and y[index+1] > half_max:
		roots.append(x[index]+(half_max-y[index])/(y[index+1]-y[index])*(x[index+1]-x[index]))
	    elif y[index] >= half_max and y[index+1] < half_max:
		roots.append(x[index]+(half_max-y[index])/(y[index+1]-y[index])*(x[index+1]-x[index]))
    elif method == 'Spline' :
	#version A / spline interpolation
	s = splrep(x, y - half_max)
	roots = sproot(s)

    if len(roots) > 2:
	x_max = 0
	for index in range(len(x)):
	    if y[index] == np.max(y):
		x_max = x[index]
	if x_max == 0:
	    return 0, 0
	else:
	    try:
	    	left = max(roots[roots[:]< x_max])
	    except:
		left = 0
	    try:
	    	right = min(roots[roots[:]> x_max])
	    except:
		right = 0
	    return (right+left)/2., abs(right - left)
    elif len(roots) < 2:
	return 0, 0
    else:
        return (roots[1] + roots[0])/2., abs(roots[1] - roots[0])


def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def linear(x, *p):
    a, b = p
    return a*x + b

def doublegauss(x, *p):
    A1, A2, mu1, mu2, sigma1, sigma2 = p
    return A1*np.exp(-(x-mu1)**2/(2.*sigma1**2)) + A2*np.exp(-(x-mu2)**2/(2.*sigma2**2))

def add_MetProjections(config, inputDataPlot, dictPlot):
    variable = 'recoMetOnGenMetProjectionPar'
    if 'LongZCorrectedRecoil_MET' in dictPlot and 'LongZCorrectedRecoil_METPhi' in dictPlot and 'genMet_Pt' in dictPlot and 'genMet_Phi' in dictPlot:
    	if not 'recoMetOnGenMetProjectionPar' in dictPlot:
            dictPlot['recoMetOnGenMetProjectionPar'] = inputDataPlot.shape[1]
            inputDataPlot = np.hstack((inputDataPlot, np.array(np.cos(inputDataPlot[:,dictPlot['genMet_Phi']]-inputDataPlot[:,dictPlot['LongZCorrectedRecoil_METPhi']])*(inputDataPlot[:,dictPlot['LongZCorrectedRecoil_MET']])).reshape(inputDataPlot.shape[0],1)))
    	if not 'recoMetOnGenMetProjectionPerp' in dictPlot:
            dictPlot['recoMetOnGenMetProjectionPerp'] = inputDataPlot.shape[1]
            inputDataPlot = np.hstack((inputDataPlot, np.array(np.sin(inputDataPlot[:,dictPlot['genMet_Phi']]-inputDataPlot[:,dictPlot['LongZCorrectedRecoil_METPhi']])*(inputDataPlot[:,dictPlot['LongZCorrectedRecoil_MET']])).reshape(inputDataPlot.shape[0],1)))


    for index in range(len(config['inputFile'])-1):
        if 'V%iLongZCorrectedRecoil_MET'%index in dictPlot and 'V%iLongZCorrectedRecoil_METPhi'%index in dictPlot and 'genMet_Pt' in dictPlot and 'genMet_Phi' in dictPlot:
            if not 'V%irecoMetOnGenMetProjectionPar'%index in dictPlot:
                dictPlot['V%irecoMetOnGenMetProjectionPar'%index] = inputDataPlot.shape[1]
                inputDataPlot = np.hstack((inputDataPlot, np.array(np.cos(inputDataPlot[:,dictPlot['genMet_Phi']]-inputDataPlot[:,dictPlot['V%iLongZCorrectedRecoil_METPhi'%index]])*(inputDataPlot[:,dictPlot['V%iLongZCorrectedRecoil_MET'%index]])).reshape(inputDataPlot.shape[0],1)))
            if not 'V%irecoMetOnGenMetProjectionPerp'%index in dictPlot:
                dictPlot['V%irecoMetOnGenMetProjectionPerp'%index] = inputDataPlot.shape[1]
                inputDataPlot = np.hstack((inputDataPlot, np.array(np.sin(inputDataPlot[:,dictPlot['genMet_Phi']]-inputDataPlot[:,dictPlot['V%iLongZCorrectedRecoil_METPhi'%index]])*(inputDataPlot[:,dictPlot['V%iLongZCorrectedRecoil_MET'%index]])).reshape(inputDataPlot.shape[0],1)))


    if 'dpfmet_Pt' in dictPlot and 'dpfmet_Phi' in dictPlot and 'genMet_Pt' in dictPlot:
    	if not 'recoPfMetOnGenMetProjectionPar' in dictPlot:
	    dictPlot['recoPfMetOnGenMetProjectionPar'] = inputDataPlot.shape[1]
            inputDataPlot = np.hstack((inputDataPlot, np.array(np.cos(inputDataPlot[:,dictPlot['dpfmet_Phi']])*(inputDataPlot[:,dictPlot['genMet_Pt']]-inputDataPlot[:,dictPlot['dpfmet_Pt']])).reshape(inputDataPlot.shape[0],1)))
    	if not 'recoPfMetOnGenMetProjectionPerp' in dictPlot:
	    dictPlot['recoPfMetOnGenMetProjectionPerp'] = inputDataPlot.shape[1]
            inputDataPlot = np.hstack((inputDataPlot, np.array(np.sin(inputDataPlot[:,dictPlot['dpfmet_Phi']])*(inputDataPlot[:,dictPlot['genMet_Pt']]-inputDataPlot[:,dictPlot['dpfmet_Pt']])).reshape(inputDataPlot.shape[0],1)))

    return inputDataPlot, dictPlot

def load_datasetcsv(config):
    #create Treevariable

    start = time.time()


    ArtusDict = {
            "mvamet" : "MVAMET_Pt",
            "mvametphi" : "MVAMET_Phi",
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
            "genMetPhi" : "genMet_Phi",
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

    if len(config['inputFile']) > 1:
	trainingheader =  ["LongZCorrectedRecoil_LongZ","LongZCorrectedRecoil_PerpZ","LongZCorrectedRecoil_Phi","PhiCorrectedRecoil_LongZ","PhiCorrectedRecoil_PerpZ","PhiCorrectedRecoil_Phi", "dmvamet_Pt", "dmvamet_Phi", "recoMetOnGenMetProjectionPar", "recoMetOnGenMetProjectionPerp","recoPfMetOnGenMetProjectionPar","recoPfMetOnGenMetProjectionPerp"]
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


    #add rojection on real Met for data from MapAnalyzer as it is not calculated there
    inputDataPlot, dictPlot = add_MetProjections(config, inputDataPlot, dictPlot)


    print(dictPlot)
    #print(inputDataPlot.shape)



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
    plt.xlabel(variablename,fontsize = 17)
    plt.ylabel('Hits',fontsize = 17)


    plt.savefig((outputdir+variablename+".png"))
    plt.clf()
    return 0


def make_ResponseCorrectedPlot(config, XRange, YStd, YResponse, bosonName, targetvariable,  minrange,maxrange, stepwidth, ptmin,ptmax, labelname = 'MVAMet', relateVar = 'p_t^Z', relateUnits = 'GeV'):

    plt.clf()
    ResCorr = YStd[:]/YResponse[:]
    plt.plot(XRange[:-1]+stepwidth/2.,ResCorr,'o')
    plt.xlabel(targetvariable,fontsize = 17)
    plt.ylabel('Resolution / Response',fontsize = 17)
    if ptmax == 0:
	  plt.savefig(config['outputDir'] + "ControlPlots/ResponseCorrected_%s_vs_%s.png" %(bosonName,targetvariable))
    else:
	  plt.savefig(config['outputDir'] + "ControlPlots/ResponseCorrected(%i<Pt<%i)_%s_vs_%s.png" %(ptmin,ptmax,bosonName,targetvariable))
    plt.figure(6)
    plt.plot(XRange[:-1]+stepwidth/2.,ResCorr,'o-',label=labelname)
    plt.figure(0)
    plt.clf()



    return

def make_ResolutionPlot(config,plotData,dictPlot, bosonName, targetvariable,  minrange=42,maxrange=0, stepwidth=0, ptmin =0,ptmax=0, labelname = 'MVAMet', relateVar = 'p_t^Z', relateUnits = 'GeV'):


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
        fitDistri = currentDistri[((currentDistri.mean()-4*currentDistri.std())<currentDistri[:]) & (currentDistri[:] <(currentDistri.mean()+4*currentDistri.std()))]

	if fitDistri.shape[0] == 0:
	    YMean[index] = 0
	    YStd[index] = 0
	else:
	    num_bins = 100

	    #YMean[index] = currentDistri.mean()
            #YStd[index] = currentDistri.std()
            #print("Left %i outlier events out for fitting"%(currentDistri.shape[0]-fitDistri.shape[0]))

            n, bins, patches = plt.hist(fitDistri, num_bins, facecolor='green', alpha=0.5)

            XCenter = (bins[:-1] + bins[1:])/2
            p0 = [1., 0., 1.]
	    try:
	            coeff, var_matrix = curve_fit(gauss,XCenter,n,p0=p0)
            except:
		    coeff = p0

            chiNormal = ((gauss(XCenter,*coeff)-n)**2/n).sum()
            chiperDOF = chiNormal/(num_bins-3)

	    FWHM_mean, FWHM_std = fwhm(XCenter,n)
            """
            coeffDouble, var_matrixDouble = curve_fit(doublegauss,XCenter,n,p0=p0)

            plt.text(meanLoc-2.5*stdLoc,0.4*(plt.ylim()[1]-plt.ylim()[0]),r'$\mu_{2Gauss}^1 = %.2f$''\n'r'$=\sigma_{2Gauss}^1= %.2f$''\n'r'$\mu_{2Gauss}^2 = %.2f$''\n'r'$=\sigma_{2Gauss}^2= %.2f$''\n'%(coeffDouble[1], coeffDouble[2], coeffDouble[4, coeffDouble[5]]),color = 'k',fontsize=16)

            print(coeffDouble)
            """
            if bosonName[0] == 'V':
		if len(config['method']) > 1:
		    if config['method'][int(bosonName[1])+1][0:5] == 'Alpha':
			YMean[index], YStd[index] = alphaFit(dictPlot, AlternativeDistri, bosonName,'Resolution', method = config['method'][int(bosonName[1])+1])
		    else:
			YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][int(bosonName[1])+1])
			YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][int(bosonName[1])+1])
		else:
		    if config['method'][0][0:5] == 'Alpha':
			YMean[index], YStd[index] = alphaFit(dictPlot, AlternativeDistri, bosonName,'Resolution', method = config['method'][0])
		    else:
			YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][0])
			YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][0])
            else:
                if config['method'][0][0:5] == 'Alpha':
		    YMean[index], YStd[index] = alphaFit(dictPlot, AlternativeDistri, bosonName,'Resolution', method = config['method'][0])
                else:
                    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][0])
                    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][0])
            meanLoc = coeff[1]
            stdLoc = abs(coeff[2])

	    alpha_mean, alpha_std = alphaFit(dictPlot, AlternativeDistri, bosonName,'Resolution')

            plt.clf()

            #y = mlab.normpdf(bins, currentDistri.mean(), currentDistri.std())
            plt.rc('font', family='serif')
            n, bins, patches = plt.hist(currentDistri, num_bins, normed=1, range=[fitDistri.mean()-3*fitDistri.std(), fitDistri.mean()+5*fitDistri.std()], facecolor='green', alpha=0.5)
            y = mlab.normpdf(bins, meanLoc, stdLoc)
            #plt.xlabel('%s at %f, mean: %f'%(targetvariable,(XRange[index+1]+XRange[index])/2,currentDistri.mean()),fontsize = 17)
            plt.xlabel(r'$U_\| - p_t^Z$ at $%s = (%i - %i)\,\mathrm{%s}$'%(relateVar,XRange[index],XRange[index+1],relateUnits),fontsize = 17)
	    plt.text(meanLoc+1.8*stdLoc,0.15*(plt.ylim()[1]-plt.ylim()[0]),r'$\mathrm{Events} = %i$''\n'r'$\mathrm{Outliers}(>4\sigma) = %.2f$%%''\n'r'$\mu_{tot} = %.3f$''\n'r'$\sigma_{tot} = %.3f$''\n'r'$\mu_{sel} = %.3f (\Delta_{tot} = %.2f\sigma_{tot}$)''\n'r'$\sigma_{sel} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\mu_{fit} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\sigma_{fit} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\chi_{pDoF}^2 = %.1f$''\n'r'$\mu_{FWHM} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\sigma_{FWHM} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\mu_{\alpha} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\sigma_{\alpha} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)'%(currentDistri.shape[0],100*(1-fitDistri.shape[0]*1./currentDistri.shape[0]),currentDistri.mean(),currentDistri.std(),fitDistri.mean(),(fitDistri.mean()-currentDistri.mean())/currentDistri.std(), fitDistri.std(),(fitDistri.std()-currentDistri.std())/currentDistri.std(), meanLoc,(meanLoc-currentDistri.mean())/currentDistri.std(),stdLoc,(stdLoc-currentDistri.std())/currentDistri.std(),chiperDOF,FWHM_mean,(FWHM_mean-currentDistri.mean())/currentDistri.std(),FWHM_std/2.355,(FWHM_std/2.355-currentDistri.std())/currentDistri.std(),alpha_mean,(alpha_mean-currentDistri.mean())/currentDistri.std(),alpha_std,(alpha_std-currentDistri.std())/currentDistri.std()),color = 'k',fontsize=16)
            #plt.title('DMean: %f , DStd: %f'%(currentDistri.mean()-coeff[1],currentDistri.std()-coeff[2]), fontsize = 20)

            #plt.ylabel('(MET Boson PT_Long)/(True Boson Pt)',fontsize = 17)
            plt.ylabel('frequency distribution',fontsize = 17)

            if ptmax == 0:
                plt.title('Resolution %s'%labelname, fontsize = 20)
                plt.plot(bins, y, 'r--')
                foldername = 'Resolution_%s_vs_%s' %(bosonName,targetvariable)
                if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
                    os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
		    bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
		    os.system(bashCommand)
                plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.png' %(foldername,foldername, index)))
		if config['method'][0][0:5] == 'Alpha':
		    alphaFit(dictPlot, AlternativeDistri, bosonName,'Resolution', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.png' %(foldername,foldername, index)), method = config['method'][0])
		else:
		    alphaFit(dictPlot, AlternativeDistri, bosonName,'Resolution', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.png' %(foldername,foldername, index)))
            else:
                plt.title('Resolution %s 'r'$(%i\,\mathrm{GeV}<p_t^Z<%i\,\mathrm{GeV})$'%(labelname,ptmin,ptmax), fontsize = 20)
                plt.plot(bins, y, 'r--')
                foldername = 'Resolution(%i<Pt<%i)_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
                Bashfoldername = 'Resolution\(%i\<Pt\<%i\)_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
                if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
                    os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
		    bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],Bashfoldername)
		    os.system(bashCommand)
                plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.png' %(foldername,foldername, index)))
		if config['method'][0][0:5] == 'Alpha':
		    alphaFit(dictPlot, AlternativeDistri, bosonName,'Resolution', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.png' %(foldername,foldername, index)), method = config['method'][0])
		else:
		    alphaFit(dictPlot, AlternativeDistri, bosonName,'Resolution', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.png' %(foldername,foldername, index)))


    plt.clf()
    plt.plot(XRange[:-1]+stepwidth/2.,YStd[:],'o-')
    plt.ylabel('(MET Boson PT_Long) - (True Boson Pt)',fontsize = 17)
    plt.xlabel(targetvariable,fontsize = 17)
    if ptmax == 0:
	plt.savefig(config['outputDir'] + "ControlPlots/Resolution_%s_vs_%s.png" %(bosonName,targetvariable))
    else:
	plt.savefig(config['outputDir'] + "ControlPlots/Resolution(%i<Pt<%i)_%s_vs_%s.png" %(ptmin,ptmax,bosonName,targetvariable))
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

	if bosonName[0] == 'V':
	    if len(config['method']) > 1:
		if config['method'][int(bosonName[1])+1][0:5] == 'Alpha':
		    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', 'Trunc')
                    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', 'Trunc')
		else:
		    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][int(bosonName[1])+1])
		    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][int(bosonName[1])+1])
	    else:
		if config['method'][0][0:5] == 'Alpha':
		    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', 'Trunc')
                    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', 'Trunc')
		else:
		    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][0])
		    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][0])
	else:
	    if config['method'][0][0:5] == 'Alpha':
		YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', 'Trunc')
		YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', 'Trunc')
	    else:
		YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][0])
		YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][0])



	if index < 12:
	    plt.clf()
	    num_bins = 50
	    n, bins, patches = plt.hist(currentDistri, num_bins, normed=1, facecolor='green', alpha=0.5)
	    y = mlab.normpdf(bins, currentDistri.mean(), currentDistri.std())
	    plt.xlabel('%s at %f, mean: %f'%(targetvariable,(XRange[index+1]+XRange[index])/2,currentDistri.mean()),fontsize = 17)
	    plt.ylabel('(MET) - (gen MET)',fontsize = 17)
	    plt.plot(bins, y, 'r--')
	    if ptmax == 0:
		foldername = 'METResolution_%s_vs_%s' %(bosonName,targetvariable)
		if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
		    os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
		    bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
		    os.system(bashCommand)
		plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.png' %(foldername,foldername, index)))
	    else:
		foldername = 'METResolution(%i<Pt<%i)_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
		Bashfoldername = 'METResolution\(%i\<Pt\<%i\)_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
		if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
		    os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
		    bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],Bashfoldername)
		    os.system(bashCommand)
		plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.png' %(foldername,foldername, index)))

    plt.clf()
    plt.plot(XRange[:-1]+stepwidth/2.,YStd[:],'o')
    plt.ylabel('(MET) - (gen MET)',fontsize = 17)
    plt.xlabel(targetvariable,fontsize = 17)
    if ptmax == 0:
	plt.savefig(config['outputDir'] + "ControlPlots/METResolution_%s_vs_%s.png" %(bosonName,targetvariable))
    else:
	plt.savefig(config['outputDir'] + "ControlPlots/METResolution(%i<Pt<%i)_%s_vs_%s.png" %(ptmin,ptmax,bosonName,targetvariable))
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
	#YMean[index] = currentDistri.mean()
	#YStd[index] = currentDistri.std()


        #YMean[index] = currentDistri.mean()
        #YStd[index] = currentDistri.std()
        fitDistri = currentDistri[((currentDistri.mean()-4*currentDistri.std())<currentDistri[:]) & (currentDistri[:] <(currentDistri.mean()+4*currentDistri.std()))]
        #print("Left %i outlier events out for fitting"%(currentDistri.shape[0]-fitDistri.shape[0]))
        if fitDistri.shape[0] == 0:
	    YMean[index] = 0
	    YStd[index] = 0
	else:
	    num_bins = 50
	    n, bins, patches = plt.hist(fitDistri, num_bins, facecolor='green', alpha=0.5)

	    XCenter = (bins[:-1] + bins[1:])/2
	    p0 = [1., 0., 1.]
	    try:
		    coeff, var_matrix = curve_fit(gauss,XCenter,n,p0=p0)
	    except:
		    coeff = p0
	    p0Double = [1., 1., 1., 0., 1., 1.]
	    #coeffDouble, var_matrixDouble = curve_fit(doublegauss,XCenter,n,p0=p0Double)

	    #print(gauss(XCenter,*coeff)-n)
	    chiNormal = ((gauss(XCenter,*coeff)-n)**2/n).sum()
	    chiperDOF = chiNormal/(num_bins-3)
	    #chiDouble = ((doublegauss(XCenter,*coeffDouble)-n)**2).sum()

	    #print("Normal Gauss: %f"%chiNormal)
	    #print("Normal Gauss: %f"%chiDouble)


	    #print(coeffDouble)


            if bosonName[0] == 'V':
                if len(config['method']) > 1:
                    if config['method'][int(bosonName[1])+1][0:5] == 'Alpha':
			YMean[index], YStd[index] = alphaFit(dictPlot, AlternativeDistri, bosonName,'Response', method = config['method'][int(bosonName[1])+1])
                    else:
                        YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][int(bosonName[1])+1])
                        YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][int(bosonName[1])+1])
                else:
                    if config['method'][0][0:5] == 'Alpha':
			YMean[index], YStd[index] = alphaFit(dictPlot, AlternativeDistri, bosonName,'Response', method = config['method'][0])
                    else:
                        YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][0])
                        YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][0])
            else:
                if config['method'][0][0:5] == 'Alpha':
		    YMean[index], YStd[index] = alphaFit(dictPlot, AlternativeDistri, bosonName,'Response', method = config['method'][0])
                else:
                    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][0])
                    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][0])

	    meanLoc = coeff[1]
	    stdLoc = abs(coeff[2])

	    FWHM_mean, FWHM_std = fwhm(XCenter,n)

	    #YMean[index] = currentDistri.mean()
	    #YStd[index] = currentDistri.std()
	    plt.clf()

	    alpha_mean, alpha_std = alphaFit(dictPlot, AlternativeDistri, bosonName,'Response')
	    #y = mlab.normpdf(bins, currentDistri.mean(), currentDistri.std())
	    plt.rc('font', family='serif')
	    n, bins, patches = plt.hist(currentDistri, num_bins, range=[fitDistri.mean()-3*fitDistri.std(), fitDistri.mean()+5*fitDistri.std()], facecolor='green', alpha=0.5)
	    #y = mlab.normpdf(bins, meanLoc, stdLoc)
	    #print("y ",y.shape)
	    y = gauss(bins,*coeff)


	    #plt.xlabel('%s at %f, mean: %f'%(targetvariable,(XRange[index+1]+XRange[index])/2,currentDistri.mean()),fontsize = 17)
	    plt.xlabel(r'$U_\| / p_t^Z$ at $%s = (%i - %i)\,\mathrm{%s}$'%(relateVar,XRange[index],XRange[index+1],relateUnits),fontsize = 17)
	    plt.text(meanLoc+1.8*stdLoc,0.15*(plt.ylim()[1]-plt.ylim()[0]),r'$\mathrm{Events} = %i$''\n'r'$\mathrm{Outliers}(>4\sigma) = %.2f$%%''\n'r'$\mu_{tot} = %.3f$''\n'r'$\sigma_{tot} = %.3f$''\n'r'$\mu_{sel} = %.3f (\Delta_{tot} = %.2f\sigma_{tot}$)''\n'r'$\sigma_{sel} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\mu_{fit} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\sigma_{fit} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\chi_{pDoF}^2 = %.1f$''\n'r'$\mu_{FWHM} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\sigma_{FWHM} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\mu_{\alpha} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\sigma_{\alpha} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)'%(currentDistri.shape[0],100*(1-fitDistri.shape[0]*1./currentDistri.shape[0]),currentDistri.mean(),currentDistri.std(),fitDistri.mean(),(fitDistri.mean()-currentDistri.mean())/currentDistri.std(), fitDistri.std(),(fitDistri.std()-currentDistri.std())/currentDistri.std(), meanLoc,(meanLoc-currentDistri.mean())/currentDistri.std(),stdLoc,(stdLoc-currentDistri.std())/currentDistri.std(),chiperDOF,FWHM_mean,(FWHM_mean-currentDistri.mean())/currentDistri.std(),FWHM_std/2.355,(FWHM_std/2.355-currentDistri.std())/currentDistri.std(),alpha_mean,(alpha_mean-currentDistri.mean())/currentDistri.std(),alpha_std,(alpha_std-currentDistri.std())/currentDistri.std()),color = 'k',fontsize=16)
	    #plt.title('DMean: %f , DStd: %f'%(currentDistri.mean()-coeff[1],currentDistri.std()-coeff[2]))

	    #plt.ylabel('(MET Boson PT_Long)/(True Boson Pt)',fontsize = 17)
	    plt.ylabel('frequency distribution',fontsize = 17)

	    if ptmax == 0:
		plt.title('Response %s'%labelname, fontsize = 20)
		plt.plot(bins, y, 'r--')
		foldername = 'Response_%s_vs_%s' %(bosonName,targetvariable)
		if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
		    os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
		    bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
		    os.system(bashCommand)
		plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.png' %(foldername,foldername, index)))
		if config['method'][0][0:5] == 'Alpha':
		    alphaFit(dictPlot, AlternativeDistri, bosonName,'Response', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.png' %(foldername,foldername, index)), method = config['method'][0])
		else:
		    alphaFit(dictPlot, AlternativeDistri, bosonName,'Response', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.png' %(foldername,foldername, index)))
	    else:
		plt.title('Response %s 'r'$(%i\,\mathrm{GeV}<p_t^Z<%i\,\mathrm{GeV})$'%(labelname,ptmin,ptmax), fontsize = 20)
		plt.plot(bins, y, 'r--')
		foldername = 'Response(%i<Pt<%i)_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
		Bashfoldername = 'Response\(%i\<Pt\<%i\)_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
		if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
		    os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
		    bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],Bashfoldername)
		    os.system(bashCommand)
		plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.png' %(foldername,foldername, index)))
		if config['method'][0][0:5] == 'Alpha':
		    alphaFit(dictPlot, AlternativeDistri, bosonName,'Response', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.png' %(foldername,foldername, index)), method = config['method'][0])
		else:
		    alphaFit(dictPlot, AlternativeDistri, bosonName,'Response', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.png' %(foldername,foldername, index)))
    plt.clf()
    plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o-')

    plt.xlabel(targetvariable,fontsize = 17)
    if ptmax == 0:
	plt.savefig(config['outputDir'] + "ControlPlots/Response_%s_vs_%s.png" %(bosonName,targetvariable))
    else:
	plt.savefig(config['outputDir'] + "ControlPlots/Response(%i<Pt<%i)_%s_vs_%s.png" %(ptmin,ptmax,bosonName,targetvariable))
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

	if bosonName[0] == 'V':
            if len(config['method']) > 1:
                if config['method'][int(bosonName[1])+1][0:5] == 'Alpha':
                    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', 'Trunc')
                    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', 'Trunc')
                else:
                    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][int(bosonName[1])+1])
                    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][int(bosonName[1])+1])
            else:
                if config['method'][0][0:5] == 'Alpha':
                    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', 'Trunc')
                    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', 'Trunc')
                else:
                    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][0])
                    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][0])
        else:
            if config['method'][0][0:5] == 'Alpha':
                YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', 'Trunc')
                YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', 'Trunc')
            else:
                YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][0])
                YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][0])

	if index < 12:
	    plt.clf()
	    num_bins = 50
	    n, bins, patches = plt.hist(currentDistri, num_bins, normed=1, facecolor='green', alpha=0.5)
	    y = mlab.normpdf(bins, currentDistri.mean(), currentDistri.std())
	    plt.xlabel('%s at %f, mean: %f'%(targetvariable,(XRange[index+1]+XRange[index])/2,currentDistri.mean()),fontsize = 17)
	    plt.ylabel('(MET)/(gen MET)',fontsize = 17)
	    plt.plot(bins, y, 'r--')
	    if ptmax == 0:
		foldername = 'METResponse_%s_vs_%s' %(bosonName,targetvariable)
		if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
		    os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
		    bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
		    os.system(bashCommand)
		plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.png' %(foldername,foldername, index)))
	    else:
		foldername = 'METResponse(%i<Pt<%i)_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
		Bashfoldername = 'METResponse\(%i\<Pt\<%i\)_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
		if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
		    os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
		    bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],Bashfoldername)
		    os.system(bashCommand)
		plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.png' %(foldername,foldername, index)))

    plt.clf()
    plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o')

    plt.xlabel(targetvariable,fontsize = 17)
    if ptmax == 0:
	plt.savefig(config['outputDir'] + "ControlPlots/METResponse_%s_vs_%s.png" %(bosonName,targetvariable))
    else:
	plt.savefig(config['outputDir'] + "ControlPlots/METResponse(%i<Pt<%i)_%s_vs_%s.png" %(ptmin,ptmax,bosonName,targetvariable))
    plt.clf()
    plt.figure(9)
    plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o-',label=labelname)

    plt.figure(0)


    return 0



def make_ResolutionPerpPlot(config,plotData,dictPlot, bosonName, targetvariable,  minrange=42,maxrange=0, stepwidth=0, ptmin =0,ptmax=0, labelname = 'MVAMet', relateVar = 'p_t^Z', relateUnits = 'GeV'):


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


	if bosonName[0] == 'V':
            if len(config['method']) > 1:
                if config['method'][int(bosonName[1])+1][0:5] == 'Alpha':
                    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', 'Trunc')
                    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', 'Trunc')
                else:
                    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][int(bosonName[1])+1])
                    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][int(bosonName[1])+1])
            else:
                if config['method'][0][0:5]== 'Alpha':
                    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', 'Trunc')
                    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', 'Trunc')
                else:
                    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][0])
                    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][0])
        else:
            if config['method'][0][0:5]  == 'Alpha':
                YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', 'Trunc')
                YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', 'Trunc')
            else:
                YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][0])
                YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][0])

	if index < 12:
	    plt.clf()
	    num_bins = 50
	    n, bins, patches = plt.hist(currentDistri, num_bins, normed=1, facecolor='green', alpha=0.5)
	    y = mlab.normpdf(bins, currentDistri.mean(), currentDistri.std())
	    plt.xlabel('%s at %f, mean: %f'%(targetvariable,(XRange[index+1]+XRange[index])/2,currentDistri.mean()),fontsize = 17)
	    plt.ylabel('MET Boson PT_Perp',fontsize = 17)
	    plt.plot(bins, y, 'r--')
	    if ptmax == 0:
		foldername = 'ResolutionPerp_%s_vs_%s' %(bosonName,targetvariable)
		if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
		    os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
		    bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
		    os.system(bashCommand)
		plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.png' %(foldername,foldername, index)))
	    else:
		foldername = 'ResolutionPerp(%i<Pt<%i)_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
		Bashfoldername = 'ResolutionPerp\(%i\<Pt\<%i\)_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
		if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
		    os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
		    bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],Bashfoldername)
		    os.system(bashCommand)
		plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.png' %(foldername,foldername, index)))

    plt.clf()
    plt.plot(XRange[:-1]+stepwidth/2.,YStd[:],'o')
    plt.ylabel('MET Boson PT_Perp',fontsize = 17)
    plt.xlabel(targetvariable,fontsize = 17)
    if ptmax == 0:
	plt.savefig(config['outputDir'] + "ControlPlots/ResolutionPerp_%s_vs_%s.png" %(bosonName,targetvariable))
    else:
	plt.savefig(config['outputDir'] + "ControlPlots/ResolutionPerp(%i<Pt<%i)_%s_vs_%s.png" %(ptmin,ptmax,bosonName,targetvariable))

    plt.figure(8)
    plt.plot(XRange[:-1]+stepwidth/2.,YStd[:],'o-',label=labelname)
    plt.figure(0)
    plt.clf()


    return

def make_ResponsePerpPlot(config, plotData,dictPlot, bosonName, targetvariable,  minrange=42,maxrange=0, stepwidth=0, ptmin=0, ptmax=0, labelname = 'MVAMet', relateVar = 'p_t^Z', relateUnits = 'GeV'):

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

	if bosonName[0] == 'V':
            if len(config['method']) > 1:
                if config['method'][int(bosonName[1])+1][0:5] == 'Alpha':
                    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', 'Trunc')
                    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', 'Trunc')
                else:
                    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][int(bosonName[1])+1])
                    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][int(bosonName[1])+1])
            else:
                if config['method'][0][0:5] == 'Alpha':
                    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', 'Trunc')
                    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', 'Trunc')
                else:
                    YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][0])
                    YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][0])
        else:
            if config['method'][0][0:5] == 'Alpha':
                YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', 'Trunc')
                YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', 'Trunc')
            else:
                YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][0])
                YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][0])

	if index < 12:
	    plt.clf()
	    num_bins = 50
	    n, bins, patches = plt.hist(currentDistri, num_bins, normed=1, facecolor='green', alpha=0.5)
	    y = mlab.normpdf(bins, currentDistri.mean(), currentDistri.std())
	    plt.xlabel('%s at %f, mean: %f'%(targetvariable,(XRange[index+1]+XRange[index])/2,currentDistri.mean()),fontsize = 17)
	    plt.ylabel('MET Boson PT_Perp',fontsize = 17)
	    plt.plot(bins, y, 'r--')
	    if ptmax == 0:
		foldername = 'ResponsePerp_%s_vs_%s' %(bosonName,targetvariable)
		if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
		    os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
		    bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
		    os.system(bashCommand)
		plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.png' %(foldername,foldername, index)))
	    else:
		foldername = 'ResponsePerp(%i<Pt<%i)_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
		Bashfoldername = 'ResponsePerp\(%i\<Pt\<%i\)_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
		if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
		    os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
		    bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],Bashfoldername)
		    os.system(bashCommand)
		plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.png' %(foldername,foldername, index)))
    plt.clf()
    plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o')

    plt.xlabel(targetvariable,fontsize = 17)
    if ptmax == 0:
	plt.savefig(config['outputDir'] + "ControlPlots/ResponsePerp_%s_vs_%s.png" %(bosonName,targetvariable))
    else:
	plt.savefig(config['outputDir'] + "ControlPlots/ResponsePerp(%i<Pt<%i)_%s_vs_%s.png" %(ptmin,ptmax,bosonName,targetvariable))
    plt.clf()
    plt.figure(7)
    plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o-',label=labelname)
    plt.figure(0)



    return



def make_ControlPlots(config, plotData,dictPlot, bosonName, targetvariable,  minrange=42,maxrange=0, stepwidth=0, ptmin=0,ptmax=0, labelname = 'MVAMet', relateVar = 'p_t^Z', relateUnits = 'GeV'):

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

    if bosonName == "LongZCorrectedRecoil" and "recoMetOnGenMetProjectionPar" in dictPlot and not plotData[:,dictPlot["recoMetOnGenMetProjectionPar"]].mean() == -999:
        make_METResponsePlot(config, plotData, dictPlot, "recoMetOnGenMetProjectionPar", targetvariable,  minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits)
        make_METResolutionPlot(config, plotData, dictPlot, "recoMetOnGenMetProjectionPar", targetvariable,  minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits)

    if bosonName[0] == "V" and "V%srecoMetOnGenMetProjectionPar"%bosonName[1] in dictPlot and not plotData[:,dictPlot["V%srecoMetOnGenMetProjectionPar"%bosonName[1]]].mean() == -999:
        make_METResponsePlot(config, plotData, dictPlot, "V%srecoMetOnGenMetProjectionPar"%bosonName[1], targetvariable,  minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits)
        make_METResolutionPlot(config, plotData, dictPlot, "V%srecoMetOnGenMetProjectionPar"%bosonName[1], targetvariable,  minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits)

    if bosonName == "recoilslimmedMETs" and "recoPfMetOnGenMetProjectionPar" in dictPlot and not plotData[:,dictPlot["recoPfMetOnGenMetProjectionPar"]].mean() == -999:
        make_METResponsePlot(config, plotData, dictPlot, "recoPfMetOnGenMetProjectionPar", targetvariable,  minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits)
        make_METResolutionPlot(config, plotData, dictPlot, "recoPfMetOnGenMetProjectionPar", targetvariable,  minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits)



    return

def make_JetStudyPlots(config, plotData, dictPlot):
    originalOutputDir = config['outputDir']
    config['outputDir'] = config['outputDir'] + '/JetStudies'
    lowerAlphaCut = config[config['activePlotting']]['lowerAlphaCut']
    upperAlphaCut = config[config['activePlotting']]['upperAlphaCut']
    print('Whole data shape: ',plotData.shape)
    lowerData = plotData[(0.1<plotData[:,dictPlot['Jet0_Pt']])]
    lowerData = plotData[(lowerAlphaCut>(plotData[:,dictPlot['Jet1_Pt']]/plotData[:,dictPlot['Jet0_Pt']]))]
    print('Lower data shape: ',lowerData.shape)

    config['outputDir'] = config['outputDir'] + '/LowerAlpha%.1f'%lowerAlphaCut
    if 'LongZCorrectedRecoil_LongZ' in dictPlot:
		bosonNameLong = 'LongZCorrectedRecoil_LongZ'
		bosonNamePerp = 'LongZCorrectedRecoil_PerpZ'
    		if not os.path.exists((config['outputDir'] + 'ControlPlots/SingleDistributions/')):
        		os.makedirs((config['outputDir'] + 'ControlPlots/SingleDistributions/'))
		YResponse = make_ResponsePlot(config, lowerData, dictPlot, bosonNameLong, 'Boson_Pt',  10,200,10,0,0,'PFMet','p_t^Z','GeV')
		XRange, YVariance = make_ResolutionPlot(config, lowerData, dictPlot, bosonNameLong, 'Boson_Pt',  10,200,10,0,0,'PFMet','p_t^Z','GeV')

    if 'recoilslimmedMETs_LongZ' in dictPlot:
		bosonNameLong = 'recoilslimmedMETs_LongZ'
		bosonNamePerp = 'recoilslimmedMETs_PerpZ'
    		if not os.path.exists((config['outputDir'] + 'ControlPlots/SingleDistributions/')):
        		os.makedirs((config['outputDir'] + 'ControlPlots/SingleDistributions/'))
		YResponse = make_ResponsePlot(config, lowerData, dictPlot, bosonNameLong, 'Boson_Pt',  10,200,10,0,0,'PFMet','p_t^Z','GeV')
		XRange, YVariance = make_ResolutionPlot(config, lowerData, dictPlot, bosonNameLong, 'Boson_Pt',  10,200,10,0,0,'PFMet','p_t^Z','GeV')


    upperData = plotData[(0.1<plotData[:,dictPlot['Jet0_Pt']])]
    upperData = plotData[(upperAlphaCut<(plotData[:,dictPlot['Jet1_Pt']]/plotData[:,dictPlot['Jet0_Pt']]))]
    print('Whole data shape: ',plotData.shape)
    print('Upper data shape: ',upperData.shape)
    config['outputDir'] = originalOutputDir + '/JetStudies'
    config['outputDir'] = config['outputDir'] + '/UpperAlpha%.1f'%upperAlphaCut
    if 'LongZCorrectedRecoil_LongZ' in dictPlot:
		bosonNameLong = 'LongZCorrectedRecoil_LongZ'
		bosonNamePerp = 'LongZCorrectedRecoil_PerpZ'
    		if not os.path.exists((config['outputDir'] + 'ControlPlots/SingleDistributions/')):
        		os.makedirs((config['outputDir'] + 'ControlPlots/SingleDistributions/'))
		YResponse = make_ResponsePlot(config, upperData, dictPlot, bosonNameLong, 'Boson_Pt',  10,200,10,0,0,'PFMet','p_t^Z','GeV')
		XRange, YVariance = make_ResolutionPlot(config, upperData, dictPlot, bosonNameLong, 'Boson_Pt',  10,200,10,0,0,'PFMet','p_t^Z','GeV')

    if 'recoilslimmedMETs_LongZ' in dictPlot:
		bosonNameLong = 'recoilslimmedMETs_LongZ'
		bosonNamePerp = 'recoilslimmedMETs_PerpZ'
    		if not os.path.exists((config['outputDir'] + 'ControlPlots/SingleDistributions/')):
        		os.makedirs((config['outputDir'] + 'ControlPlots/SingleDistributions/'))
		YResponse = make_ResponsePlot(config, upperData, dictPlot, bosonNameLong, 'Boson_Pt',  10,200,10,0,0,'PFMet','p_t^Z','GeV')
		XRange, YVariance = make_ResolutionPlot(config, upperData, dictPlot, bosonNameLong, 'Boson_Pt',  10,200,10,0,0,'PFMet','p_t^Z','GeV')

def make_MoreBDTPlots(config, plotData, dictPlot):

    num_bins = 50

    if not os.path.exists((config['outputDir'] + 'CustomPlots/')):
        os.makedirs((config['outputDir'] + 'CustomPlots/'))


    if 'recoilslimmedMETs_Phi' in dictPlot and 'Boson_Phi' in dictPlot and 'PhiCorrectedRecoil_Phi' in dictPlot:
        DPhiPFBoson = plotData[:,dictPlot["Boson_Phi"]]-plotData[:,dictPlot["recoilslimmedMETs_Phi"]]+np.pi - 2.*np.pi*((plotData[:,dictPlot["Boson_Phi"]]-plotData[:,dictPlot["recoilslimmedMETs_Phi"]])>0)
        DPhiMVABoson = plotData[:,dictPlot["Boson_Phi"]]-plotData[:,dictPlot["PhiCorrectedRecoil_Phi"]]+np.pi - 2.*np.pi*((plotData[:,dictPlot["Boson_Phi"]]-plotData[:,dictPlot["PhiCorrectedRecoil_Phi"]])>0)

        plt.clf()
        n, bins, patches = plt.hist(DPhiPFBoson, num_bins, facecolor='green', label=(r'$\Delta(\phi_{PF},\phi_{Z})$'))
        plt.legend(loc='best', numpoints=1)
    	plt.xlabel(r'$\phi$',fontsize = 17)
    	plt.ylabel('Frequency distribution',fontsize = 17)
        plt.savefig(config['outputDir'] + "/CustomPlots/DPhiPFBoson.png")

        plt.clf()
        n, bins, patches = plt.hist(DPhiMVABoson, num_bins, facecolor='green', label=(r'$\Delta(\phi_{MVA},\phi_{Z})$'))

        plt.legend(loc='best', numpoints=1)
    	plt.xlabel(r'$\phi$',fontsize = 17)
    	plt.ylabel('Frequency distribution',fontsize = 17)
        plt.savefig(config['outputDir'] + "/CustomPlots/DPhiMVABoson.png")

        plt.clf()

        n, bins, patches = plt.hist([DPhiMVABoson,DPhiPFBoson], num_bins, label=[r'$\Delta(\phi_{MVA},\phi_{Z})$',r'$\Delta(\phi_{PF},\phi_{Z})$'])

        plt.legend(loc='best', numpoints=1)
        plt.xlabel(r'$\phi$',fontsize = 17)
        plt.ylabel('Frequency distribution',fontsize = 17)
        plt.savefig(config['outputDir'] + "/CustomPlots/DPhiBothBoson.png")

        plt.clf()


    if 'PhiCorrectedRecoil_LongZ' in dictPlot and 'Boson_Pt' in dictPlot:
        lowerPtCut = 150
        upperPtCut = 160

        plotDataCut = plotData[(lowerPtCut<plotData[:,dictPlot['Boson_Pt']]) & (upperPtCut>plotData[:,dictPlot['Boson_Pt']])]
        BosonPtDistri = plotDataCut[:,dictPlot["Boson_Pt"]]
        PhiPtDistri = -plotDataCut[:,dictPlot["PhiCorrectedRecoil_LongZ"]]
        TargetDistri = BosonPtDistri/PhiPtDistri

        plt.clf()
        n, bins, patches = plt.hist(BosonPtDistri, num_bins, facecolor='green', label=(r'$p_t^Z$'))
        plt.legend(loc='best', numpoints=1)
        plt.xlabel(r'$p_t$ in GeV',fontsize = 17)
        plt.ylabel('Frequency distribution',fontsize = 17)
        plt.savefig(config['outputDir'] + "/CustomPlots/BosonPt_%i_to_%i.png"%(lowerPtCut,upperPtCut))

        plt.clf()
        n, bins, patches = plt.hist(PhiPtDistri, num_bins, facecolor='green', label=(r'$\mathrm{MVAMet}_{\Phi}$'), range=[lowerPtCut-45,upperPtCut+45])
        plt.legend(loc='best', numpoints=1)
        plt.xlabel(r'$p_t$ in GeV',fontsize = 17)
        plt.ylabel('Frequency distribution',fontsize = 17)
        plt.savefig(config['outputDir'] + "/CustomPlots/PhiCorrectedPt_%i_to_%i.png"%(lowerPtCut,upperPtCut))

        plt.clf()
        n, bins, patches = plt.hist([BosonPtDistri,PhiPtDistri], num_bins, label=[r'$p_t^Z$',r'$\mathrm{MVAMet}_{\Phi}$'], range=[lowerPtCut-45,upperPtCut+45])

        plt.legend(loc='best', numpoints=1)
        plt.xlabel(r'$p_t$ in GeV',fontsize = 17)
        plt.ylabel('Frequency distribution',fontsize = 17)
        plt.savefig(config['outputDir'] + "/CustomPlots/BosonAndPhiPt_%i_to%i.png"%(lowerPtCut,upperPtCut))

        plt.clf()
        n, bins, patches = plt.hist(TargetDistri, num_bins, facecolor='green', label=(r'Target = $p_t^Z$ / $\mathrm{MVAMet}_{\Phi}$'), range=[0.5,2])
        plt.legend(loc='best', numpoints=1)
        plt.xlabel(r'Target',fontsize = 17)
        plt.ylabel('Frequency distribution',fontsize = 17)
        plt.text(0.70*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.70*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$\mu_{empiric} = %.3f$''\n'%TargetDistri.mean(),color = 'k',fontsize=16)
        plt.savefig(config['outputDir'] + "/CustomPlots/Target_%i_to_%i.png"%(lowerPtCut,upperPtCut))

        plt.clf()


    if 'recoilslimmedMETs_LongZ' in dictPlot and 'Boson_Pt' in dictPlot and 'LongZCorrectedRecoil_LongZ' in dictPlot and 'LongZCorrectedRecoil_PerpZ' in dictPlot:
        plt.clf()
        MVARecoilDiffDistri = np.sqrt((plotData[:,dictPlot["LongZCorrectedRecoil_LongZ"]]+plotData[:,dictPlot["Boson_Pt"]])**2+(plotData[:,dictPlot["LongZCorrectedRecoil_PerpZ"]])**2)
        PFRecoilDiffDistri = np.sqrt((plotData[:,dictPlot["recoilslimmedMETs_LongZ"]]+plotData[:,dictPlot["Boson_Pt"]])**2 +(plotData[:,dictPlot["recoilslimmedMETs_PerpZ"]])**2)
        n, bins, patches = plt.hist([PFRecoilDiffDistri,MVARecoilDiffDistri], num_bins, label=[r'$Reco_{Recoil}$',r'$MVA_{Recoil}$'], range=[0,80])
        plt.xlabel(r"MET in GeV",fontsize = 17)
        plt.ylabel("Frequency Distribution", fontsize = 17)
        plt.legend(loc='best', numpoints=1, fontsize = 17)
        plt.text(0.70*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.50*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$\mu_{Reco} = %.3f$''\n'r'$\mu_{MVA} = %.3f$''\n'%(PFRecoilDiffDistri.mean(),MVARecoilDiffDistri.mean()),color = 'k',fontsize=17)
        plt.savefig(config['outputDir'] + "/CustomPlots/Diff_Recoil_BosonPt.png")
        plt.clf()

    print("Custom plots created")



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
    plt.xlabel('Variance %s from true Boson Phi (%i < Boson Pt < %i)GeV. MSE: %f'%(xlabelname, ptmin, ptmax, MSE),fontsize = 17)
    plt.ylabel('Entries',fontsize = 17)
    plt.savefig(config['outputDir'] + "/CustomPlots/PhiVariance%s_%ito%iPt.png"%(xlabelname,ptmin,ptmax))
    plt.clf()
    print('MSE %s (%i<Pt<%i): '%(xlabelname,ptmin,ptmax),MSE)

    # normal distribution center at x=0 and y=5
    plt.hist2d(plotData[:,dictPlot['Boson_Phi']], histDataPhi,bins = 80, norm=LogNorm())
    #plt.ylim([-0.25,0.25])
    plt.xlabel('Boson Phi (%i < Boson Pt < %i)GeV'%(ptmin, ptmax),fontsize = 17)
    plt.ylabel('Variance of (Prediction-Target) %s'%xlabelname,fontsize = 17)
    plt.colorbar()
    plt.savefig(config['outputDir'] + "/CustomPlots/Variance2D_%s(%i<Pt<%i).png"%(xlabelname,ptmin,ptmax))
    plt.clf()




def plot_results(config, plotData, dictPlot):

    #plotData = plotData[0==plotData[:,dictPlot['NCleanedJets']],:]








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

    if plotconfig['plotAdditionalCustomPlots']:
        make_MoreBDTPlots(config, plotData, dictPlot)
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
                make_PhiVariancePlot(config, slicedData,dictPlot,'PhiCorrectedRecoil_Phi', min, bosonmax[i], 'BDT Phi')


            #BDT

            if 'LongZCorrectedRecoil_LongZ' in dictPlot:
                if 'inputLabel' in config:
                    make_ControlPlots(config, slicedData, dictPlot, 'LongZCorrectedRecoil', 'NVertex', 5,40,5,min,bosonmax[i],config['inputLabel'][0],'\# PV','')
                else:
                    make_ControlPlots(config, slicedData, dictPlot, 'LongZCorrectedRecoil', 'NVertex', 5,40,5,min,bosonmax[i],'MVAMet','\# PV','')


            for index in range(len(config['inputFile'])-1):
                if 'V%iLongZCorrectedRecoil_LongZ'%index in dictPlot:
                    if 'inputLabel' in config:
                        make_ControlPlots(config, slicedData, dictPlot, 'V%iLongZCorrectedRecoil'%index, 'NVertex', 5,40,5,min,bosonmax[i],config['inputLabel'][index+1],'\# PV','')
                    else:
                        make_ControlPlots(config, slicedData, dictPlot, 'V%iLongZCorrectedRecoil'%index, 'NVertex', 5,40,5,min,bosonmax[i],'MVAMet %i'%(index+2),'\# PV','')


       	    if plotconfig['plotPhiCorrected']:
                if 'PhiCorrectedRecoil_LongZ' in dictPlot:
                    if 'inputLabel' in config:
                        make_ControlPlots(config, slicedData, dictPlot, 'PhiCorrectedRecoil', 'NVertex', 5,40,5,min,bosonmax[i],config['inputLabel'][0],'\# PV','')
                    else:
                        make_ControlPlots(config, slicedData, dictPlot, 'PhiCorrectedRecoil', 'NVertex', 5,40,5,min,bosonmax[i],'MVAMet','\# PV','')


            for index in range(len(config['inputFile'])-1):
                if 'V%iPhiCorrectedRecoil_LongZ'%index in dictPlot:
                    if 'inputLabel' in config:
                        make_ControlPlots(config, slicedData, dictPlot, 'V%iPhiCorrectedRecoil'%index, 'NVertex', 5,40,5,min,bosonmax[i],config['inputLabel'][index+1],'\# PV','')
                    else:
                        make_ControlPlots(config, slicedData, dictPlot, 'V%iPhiCorrectedRecoil'%index, 'NVertex', 5,40,5,min,bosonmax[i],'MVAMet %i'%(index+2),'\# PV','')


        #slimmedMet
        if 'recoilslimmedMETs_LongZ' in dictPlot:
            make_ControlPlots(config, slicedData, dictPlot, 'recoilslimmedMETs', 'NVertex', 5,40,5,min,bosonmax[i],'PFMet','\# PV','')


        #Puppi-Met
        if plotconfig['plotPuppiPerformance']:
            if 'recoilslimmedMETsPuppi_LongZ' in dictPlot:
                make_ControlPlots(config, slicedData, dictPlot, 'recoilslimmedMETsPuppi', 'NVertex', 5,40,5,min,bosonmax[i],'PuppiMet','\# PV','')


        plt.clf()
        plt.figure(4)
        plt.xlabel(r'#PV ($%i\,\mathrm{GeV} < p_t^Z < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]),fontsize = 17)
        plt.ylabel(r'$<U_\| / p_t^Z>$',fontsize = 17)
        plt.title('Response 'r'$U_\|$'' vs 'r'#PV', fontsize = 20)
        legend = plt.legend(loc='best', shadow=True, numpoints=1)
        plt.plot((plt.xlim()[0], plt.xlim()[1]), (1, 1), 'k--')
        if config['fixedRange']:
            plt.ylim([0.5,1.1])
        plt.savefig(config['outputDir'] + 'Response_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()

        plt.figure(5)
        plt.xlabel(r'#PV ($%i\,\mathrm{GeV} < p_t^Z < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]),fontsize = 17)
        plt.ylabel(r'$\sigma(<U_\| - p_t^Z>)$'' in GeV',fontsize = 17)
        plt.title('Resolution 'r'$U_\|$'' vs 'r'#PV', fontsize = 20)
        legend = plt.legend(loc='best', shadow=True, numpoints=1)
        if config['fixedRange']:
            plt.ylim([10,40])
        plt.savefig(config['outputDir'] + 'Resolution_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()

        plt.figure(6)
        plt.xlabel(r'#PV ($%i\,\mathrm{GeV} < p_t^Z < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]),fontsize = 17)
        plt.ylabel('Resolution / Response',fontsize = 17)
        plt.title('Response Corrected vs 'r'#PV', fontsize = 20)
        legend = plt.legend(loc='best', shadow=True, numpoints=1)
        if config['fixedRange']:
            plt.ylim([10,40])
        plt.savefig(config['outputDir'] + 'ResponseCorrected_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()

        plt.figure(7)
        legend = plt.legend(loc='best', shadow=True, numpoints=1)
        plt.xlabel(r'#PV ($%i\,\mathrm{GeV} < p_t^Z < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]),fontsize = 17)
        plt.ylabel(r'$<U_\bot>$',fontsize = 17)
        plt.title('Response 'r'$U_\bot$'' vs 'r'#PV', fontsize = 20)
        if config['fixedRange']:
            plt.ylim([-1,1])
        plt.savefig(config['outputDir'] + 'ResponsePerp_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()

        plt.figure(8)
        plt.xlabel(r'#PV ($%i\,\mathrm{GeV} < p_t^Z < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]),fontsize = 17)
        plt.ylabel(r'$\sigma(<U_\bot>)$',fontsize = 17)
        plt.title('Resolution 'r'$U_\bot$'' vs 'r'#PV', fontsize = 20)
        legend = plt.legend(loc='best', shadow=True, numpoints=1)
        if config['fixedRange']:
            plt.ylim([10,40])
        plt.savefig(config['outputDir'] + 'ResolutionPerp_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()

        plt.figure(9)
        plt.xlabel(r'#PV ($%i\,\mathrm{GeV} < p_t^Z < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]),fontsize = 17)
        plt.ylabel(r'$E_{t\|}^{miss}/E_{t,gen}^{miss}$',fontsize = 17)
        legend = plt.legend(loc='best', shadow=True, numpoints=1)
        if config['fixedRange']:
            plt.ylim([0,2])
        plt.savefig(config['outputDir'] + 'METResponse_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()

        plt.figure(10)
        plt.ylabel('std(MET_Long/genMET)',fontsize = 17)
        plt.ylabel(r'$\sigma(E_{t\|}^{miss}-E_{t,gen}^{miss})$',fontsize = 17)
        plt.xlabel(r'#PV ($%i\,\mathrm{GeV} < p_t^Z < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]),fontsize = 17)
        legend = plt.legend(loc='best', shadow=True, numpoints=1)
        if config['fixedRange']:
            plt.ylim([10,50])
        plt.savefig(config['outputDir'] + 'METResolution_(%i<Pt<%i)_vs_NVertex'%(min,bosonmax[i]))
        plt.clf()
        plt.figure(0)


	#additional BDTPlots
        if plotconfig['plotAdditionalBDTPlots']:
            make_PhiVariancePlot(config, slicedData,dictPlot,'recoilslimmedMETs_Phi', min, bosonmax[i], 'PF Phi')

            make_PhiVariancePlot(config, slicedData,dictPlot,'recoilslimmedMETsPuppi_Phi', min, bosonmax[i],'PUPPI Phi')



    #Boson PT

    #BDT Performance
    if plotconfig['plotBDTPerformance']:
	#MVAMet
        if 'LongZCorrectedRecoil_LongZ' in dictPlot:
            #If labeling is activated
            if 'inputLabel' in config:
                make_ControlPlots(config, plotData, dictPlot, 'LongZCorrectedRecoil', 'Boson_Pt', 10,200,10,0,0,config['inputLabel'][0],'p_t^Z','GeV')
            else:
                make_ControlPlots(config, plotData, dictPlot, 'LongZCorrectedRecoil', 'Boson_Pt', 10,200,10,0,0,'MVAMet','p_t^Z','GeV')

            #for additional inputs
            for index in range(len(config['inputFile'])-1):
                if 'V%iLongZCorrectedRecoil_LongZ'%index in dictPlot:
                    if 'inputLabel' in config:
                        make_ControlPlots(config, plotData, dictPlot, 'V%iLongZCorrectedRecoil'%index, 'Boson_Pt', 10,200,10,0,0,config['inputLabel'][index+1],'p_t^Z','GeV')
                    else:
                        make_ControlPlots(config, plotData, dictPlot, 'V%iLongZCorrectedRecoil'%index, 'Boson_Pt', 10,200,10,0,0,'MVAMet %i'%(index+2),'p_t^Z','GeV')

        #Phi Corrected MVAMET Plots
        if plotconfig['plotPhiCorrected']:
            if 'PhiCorrectedRecoil_LongZ' in dictPlot:
                if 'inputLabel' in config:
                    make_ControlPlots(config, plotData, dictPlot, 'PhiCorrectedRecoil', 'Boson_Pt', 10,200,10,0,0,config['inputLabel'][0],'p_t^Z','GeV')
                else:
                    make_ControlPlots(config, plotData, dictPlot, 'PhiCorrectedRecoil', 'Boson_Pt', 10,200,10,0,0,'MVAMet','p_t^Z','GeV')

                for index in range(len(config['inputFile'])-1):
                    if 'V%iPhiCorrectedRecoil_LongZ'%index in dictPlot:
                        if 'inputLabel' in config:
                            make_ControlPlots(config, plotData, dictPlot, 'V%iPhiCorrectedRecoil'%index, 'Boson_Pt', 10,200,10,0,0,config['inputLabel'][index+1],'p_t^Z','GeV')
                        else:
                            make_ControlPlots(config, plotData, dictPlot, 'V%iPhiCorrectedRecoil'%index, 'Boson_Pt', 10,200,10,0,0,'MVAMet %i'%(index+2),'p_t^Z','GeV')

    #PF Met Control Plots Pt
    if 'recoilslimmedMETs_LongZ' in dictPlot:
        make_ControlPlots(config, plotData, dictPlot, 'recoilslimmedMETs', 'Boson_Pt', 10,200,10,0,0,'PFMet','p_t^Z','GeV')

    #Puppi Met Control Plots Pt
    if plotconfig['plotPuppiPerformance']:
        if 'recoilslimmedMETsPuppi_LongZ' in dictPlot:
            make_ControlPlots(config, plotData, dictPlot, 'recoilslimmedMETsPuppi', 'Boson_Pt', 10,200,10,0,0,'PuppiMet','p_t^Z','GeV')




    plt.figure(4)
    legend = plt.legend(loc='best', shadow=True, numpoints=1)
    #legend.get_frame().set_alpha(0.5)
    plt.ylabel(r'$<U_\| / p_t^Z>$', fontsize = 17)
    plt.xlabel(r'$p_t^Z$'' in GeV', fontsize = 17)
    plt.title('Response 'r'$U_\|$'' vs 'r'$p_t^Z$', fontsize= 20)
    plt.plot((plt.xlim()[0], plt.xlim()[1]), (1, 1), 'k--')
    if config['fixedRange']:
        plt.ylim([0.5,1.1])
    plt.savefig(config['outputDir'] + 'Response_vs_BosonPt')
    plt.clf()
    plt.figure(5)
    plt.xlabel(r'$p_t^Z$'' in GeV', fontsize = 17)
    plt.ylabel(r'$\sigma(<U_\| - p_t^Z>)$'' in GeV', fontsize = 17)
    plt.title('Resolution 'r'$U_\|$'' vs 'r'$p_t^Z$', fontsize = 20)
    legend = plt.legend(loc='best', shadow=True, numpoints=1)
    #legend.get_frame().set_alpha(0.5)
    if config['fixedRange']:
        plt.ylim([10,40])
    plt.savefig(config['outputDir'] + 'Resolution_vs_BosonPt')
    plt.clf()
    plt.figure(6)
    plt.xlabel(r'$p_t^Z$'' in GeV', fontsize = 17)
    plt.ylabel('Resolution / Response',fontsize = 17)
    plt.title('Response Corrected vs 'r'$p_t^Z$', fontsize = 20)
    legend = plt.legend(loc='best', shadow=True, numpoints=1)
    if config['fixedRange']:
        plt.ylim([10,40])
    #legend.get_frame().set_alpha(0.5)
    plt.savefig(config['outputDir'] + 'ResponseCorrected_vs_BosonPt')
    plt.clf()
    plt.figure(7)
    legend = plt.legend(loc='best', shadow=True, numpoints=1)
    #legend.get_frame().set_alpha(0.5)
    plt.xlabel(r'$p_t^Z$'' in GeV', fontsize = 17)
    plt.ylabel(r'$<U_\bot>$',fontsize = 17)
    plt.title('Response 'r'$U_\bot$'' vs 'r'$p_t^Z$', fontsize = 20)
    if config['fixedRange']:
        plt.ylim([-1,1])
    plt.savefig(config['outputDir'] + 'ResponsePerp_vs_BosonPt')
    plt.clf()
    plt.figure(8)
    plt.xlabel(r'$p_t^Z$'' in GeV',fontsize = 17)
    plt.ylabel(r'$\sigma(<U_\bot>)$',fontsize = 17)
    plt.title('Resolution 'r'$U_\bot$'' vs 'r'$p_t^Z$', fontsize = 20)
    legend = plt.legend(loc='best', shadow=True, numpoints=1)
    #legend.get_frame().set_alpha(0.5)
    if config['fixedRange']:
        plt.ylim([10,40])
    plt.savefig(config['outputDir'] + 'ResolutionPerp_vs_BosonPt')
    plt.clf()
    plt.figure(9)
    legend = plt.legend(loc='best', shadow=True, numpoints=1)
    #legend.get_frame().set_alpha(0.5)
    plt.xlabel(r'$p_t^Z$'' in GeV',fontsize = 17)
    #plt.ylabel('MET_Long/genMET',fontsize = 17)
    plt.ylabel(r'$E_{t\|}^{miss}/E_{t,gen}^{miss}$',fontsize = 17)
    #plt.ylabel(r'$\ensuremath{{\not\mathrel{E}}_T}$',fontsize = 17)
    if config['fixedRange']:
        plt.ylim([0,2])
    plt.savefig(config['outputDir'] + 'METResponse_vs_BosonPt')
    plt.clf()
    plt.figure(10)
    plt.xlabel(r'$p_t^Z$'' in GeV',fontsize = 17)
    plt.ylabel(r'$\sigma(E_{t\|}^{miss}-E_{t,gen}^{miss})$',fontsize = 17)
    legend = plt.legend(loc='best', shadow=True, numpoints=1)
    #legend.get_frame().set_alpha(0.5)
    if config['fixedRange']:
        plt.ylim([10,50])
    plt.savefig(config['outputDir'] + 'METResolution_vs_BosonPt')
    plt.clf()
    plt.figure(0)


    #Jet study plots, can be enabled in the config
    if plotconfig['plotJetStudyPlots'] and 'Jet0_Pt' in dictPlot and 'Jet1_Pt' in dictPlot:
        make_JetStudyPlots(config, plotData, dictPlot)


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
    else:
        print("Size of dataset: %i"%plotData.shape[0])

    print(plotData.shape)

    # perform plotting
    plot_results(config, plotData, dictPlot)

    #save the config as json in the plots folder
    with open(config['outputDir'] + 'config_%s.json'%os.path.basename(config['outputDir'][:-1]), 'w') as fp:
        json.dump(config, fp, sort_keys=True, indent=4)

    #exporting of the plots to ekpwww/nzaeh
    if 'export' in config:
        bashCommand = 'cp -r %s /ekpwww/nzaeh/public_html/'%config['outputDir']
        os.system(bashCommand)

    #copying index.php which allows viewing all plots at once on ekpwww by opening php.index
    bashCommand = 'cp index.php /ekpwww/nzaeh/public_html/%s/'%os.path.basename(config['outputDir'][:-1])
    os.system(bashCommand)


    print('Plots exported to ekpwww!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make MVAMet control plots.')
    parser.add_argument('-p', '--plottingconfig', default='../configs/config.json', help='Path to configurations file')
    parser.add_argument('-i', '--inputfile',nargs='+', default='', help='[optional] Inputfile(s) from which to create the plots from')
    parser.add_argument('-l', '--inputlabel',nargs='+', default='', help='[optional] Inputlabelname(s) to use in plots')
    parser.add_argument('-o', '--outputfolder', default='', help='[optional] Foldername in which to store the plots in')
    parser.add_argument('-c', '--constraints', default='', help='[optional] Constraints to data. E.g.: 50<=plotData[:,dictPlot["Boson_Pt"]] & 50<=plotData[:,dictPlot["recoilslimmedMETs_LongZ"]]')
    parser.add_argument('-e', '--export', dest='export', action='store_true', help='[optional] Exports plots to ekpwww after creating them')
    parser.add_argument('-r', '--fixedRange', dest='fixedRange', action='store_true', help='[optional] Fix ranges to default ranges in order to facilitate comparison between plots')
    parser.add_argument('-m', '--method',nargs='+', default=['Trunc'], help='[optional] Change method(s) to obtain Values. [Empiric, Trunc, Fit, FWHM, Alpha, FWHMSpline, AlphaExcl, AlphaExclInt]')
    parser.set_defaults(export=False)
    args = parser.parse_args()
    print('Used Configfile: ',args.plottingconfig)



    with open(args.plottingconfig, 'r') as f:
        config = json.load(f)

    if args.export:
        print('Exporting files to ekpwww afterwards.')
        config['export'] = True

    if not 'method' in config:
        config['method'] = args.method

    if config['method'] == 'Alpha':
        config['method'] = 'AlphaInclInt'

    if not args.inputfile == '':
        config['inputFile'] = args.inputfile

    if args.fixedRange:
        config['fixedRange'] = True
        print('Ranges will be fixed to predefined ranges.')

    if not args.inputlabel == '':
        config['inputLabel'] = args.inputlabel

    if not args.outputfolder == '':
        config['outputDir'] = "../plots/" + args.outputfolder + "/"

    print('Saving in %s'%config['outputDir'])

    if len(args.method) == 1:
        print('Used method: %s'%config['method'][0])
        if len(args.inputfile) == 1:
            if not args.constraints == '':
                config['constraints'] = args.constraints
                print('Inputfile: %s with constraint %s'%(config['inputFile'][0],args.constraints))
            else:
                print('Inputfile: %s'%config['inputFile'][0])
        else:
            for names in config['inputFile']:
                if not args.constraints == '':
                    config['constraints'] = args.constraints
                    print('Inputfile: %s with constraint %s'%(names,args.constraints))
                else:
                    print('Inputfile: %s'%names)
    else:
        for index in range(len(config['inputFile'])):
            if not args.constraints == '':
                config['constraints'] = args.constraints
                print('Inputfile: %s with constraint %s and method %s'%(args.inputfile[index],args.constraints,config['method'][index]))
            else:
                print('Inputfile: %s with method %s'%(config['inputFile'][index],config['method'][index]))

    main(config)


