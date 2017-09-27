#!/usr/bin/env python


from __future__ import print_function

import sys
import os
import os.path
import time
#import ROOT

import re
import time
import csv
import numpy as np
import math
import argparse
import pickle

import json

import matplotlib.mlab as mlab
import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['savefig.dpi'] = 500
#mpl.rcParams['figure.dpi'] = 900
mpl.rcParams['xtick.labelsize'] = 17
mpl.rcParams['ytick.labelsize'] = 17
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
def quantityComparison(config, currentDistri, dictPlot, AlternativeDistri, bosonName):

	num_bins = 50


	empiricMean = get_Quantity(currentDistri, dictPlot, 'Mean', 'Empiric')
	empiricStd = get_Quantity(currentDistri, dictPlot, 'Std', 'Empiric')

	truncMean = get_Quantity(currentDistri, dictPlot, 'Mean', 'Trunc')
	truncStd = get_Quantity(currentDistri, dictPlot, 'Std', 'Trunc')

	fitMean = get_Quantity(currentDistri, dictPlot, 'Mean', 'Fit')
	fitStd = get_Quantity(currentDistri, dictPlot, 'Std', 'Fit')

	fwhmMean = get_Quantity(currentDistri, dictPlot, 'Mean', 'FWHM')
	fwhmStd = get_Quantity(currentDistri, dictPlot, 'Std', 'FWHM')

	quantileMean = get_Quantity(currentDistri, dictPlot, 'Mean', 'Quantile')
	quantileStd = get_Quantity(currentDistri, dictPlot, 'Std', 'Quantile')

	plt.clf()

	plt.rc('font', family='serif')
	n, bins, patches = plt.hist(currentDistri[:,0], num_bins, normed= True, range=[truncMean-3*truncStd,truncMean+5*truncStd], facecolor='green', alpha=0.5,weights=currentDistri[:,1])

	if not fitStd == 1:
		y = mlab.normpdf(bins, fitMean, fitStd)
		plt.plot(bins, y, 'r--')

	truncEvents = get_Quantity(currentDistri, dictPlot, 'Mean', 'Events')

	#if "plotAlphaMethod" in config[config['activePlotting']]:
		#if config[config['activePlotting']]["plotAlphaMethod"]:
			#alphaMean, alphaStd = alphaFit(dictPlot, AlternativeDistri, bosonName,'Resolution')
			#plt.text(truncMean+1.8*truncStd,0.15*(plt.ylim()[1]-plt.ylim()[0]),r'$\mathrm{Events} = %i$''\n'r'$\mathrm{Outliers}(>4\sigma) = %.2f$%%''\n'r'$\mu_{tot} = %.3f$''\n'r'$\sigma_{tot} = %.3f$''\n'r'$\mu_{sel} = %.3f (\Delta_{tot} = %.2f\sigma_{tot}$)''\n'r'$\sigma_{sel} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\mu_{fit} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\sigma_{fit} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\mu_{FWHM} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\sigma_{FWHM} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\mu_{\alpha} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)''\n'r'$\sigma_{\alpha} = %.3f (\Delta_{tot} = %.2f \sigma_{tot}$)'%(currentDistri.shape[0],100*(1-truncEvents*1./currentDistri.shape[0]),empiricMean,empiricStd,truncMean,(truncMean-empiricMean)/empiricStd, truncStd,(truncStd-empiricStd)/empiricStd, fitMean,(fitMean-empiricMean)/empiricStd,fitStd,(fitStd-empiricStd)/empiricStd,fwhmMean,(fwhmMean-empiricMean)/empiricStd,fwhmStd/2.355,(fwhmStd/2.355-empiricStd)/empiricStd,alphaMean,(alphaMean-empiricMean)/empiricStd,alphaStd,(alphaStd-empiricStd)/empiricStd),color = 'k',fontsize=16)
	#else:
		#print("TEST!")
	#plt.text(truncMean+1.8*truncStd,0.21*(plt.ylim()[1]-plt.ylim()[0]),r'$\mathrm{Events} = %i$''\n'r'$\mathrm{Outliers}(>4\sigma) = %.2f$%%''\n'r'$\mu_{tot} = %.3f$''\n'r'$\sigma_{tot} = %.3f$''\n'r'$\mu_{sel} = %.3f$ ''\n'r'$\sigma_{sel} = %.3f$''\n'r'$\mu_{fit} = %.3f$''\n'r'$\sigma_{fit} = %.3f $''\n'r'$\mu_{FWHM} = %.3f $''\n'r'$\sigma_{FWHM} = %.3f $''\n'r'$\mu_{Qtl} = %.3f $''\n'r'$\sigma_{Qtl} = %.3f$'%(currentDistri.shape[0],100*(1-truncEvents*1./currentDistri.shape[0]),empiricMean,empiricStd,truncMean, truncStd, fitMean,fitStd,fwhmMean,fwhmStd,quantileMean,quantileStd),color = 'k',fontsize=16)
	plt.text(truncMean+1.8*truncStd,0.21*(plt.ylim()[1]-plt.ylim()[0]),r'$\mathrm{Events} = %i$''\n'r'$\mathrm{Outliers}(>4\sigma) = %.2f$%%''\n'r'$\mu_{tot} = %.3f$''\n'r'$\mu_{sel} = %.3f$ ''\n'r'$\mu_{fit} = %.3f$''\n'r'$\mu_{FWHM} = %.3f $''\n'r'$\mu_{Qtl} = %.3f $''\n'r'$\sigma_{tot} = %.3f$''\n'r'$\sigma_{sel} = %.3f$ ''\n'r'$\sigma_{fit} = %.3f$''\n'r'$\sigma_{FWHM} = %.3f $''\n'r'$\sigma_{Qtl} = %.3f $'%(currentDistri.shape[0],100*(1-truncEvents*1./currentDistri.shape[0]),empiricMean,truncMean,fitMean,fwhmMean,quantileMean,empiricStd,truncStd,fitStd,fwhmStd,quantileStd),color = 'k',fontsize=16)
	plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
	#plt.text(truncMean+1.5*truncStd,0.21*(plt.ylim()[1]-plt.ylim()[0]),r'$\mathrm{Events} = %i$''\n'r'$\mathrm{Outliers}(>4\sigma) = %.2f$%%''\n'r'$\mu_{tot} = %.3f$''\n'r'$\sigma_{tot} = %.3f$''\n'r'$\mu_{sel} = %.3f (\Delta = %.2f\sigma_{tot}$)''\n'r'$\sigma_{sel} = %.3f (\Delta = %.2f \sigma_{tot}$)''\n'r'$\mu_{fit} = %.3f (\Delta = %.2f \sigma_{tot}$)''\n'r'$\sigma_{fit} = %.3f (\Delta = %.2f \sigma_{tot}$)''\n'r'$\mu_{FWHM} = %.3f (\Delta = %.2f \sigma_{tot}$)''\n'r'$\sigma_{FWHM} = %.3f (\Delta = %.2f \sigma_{tot}$)''\n'r'$\mu_{Qtl} = %.3f (\Delta = %.2f \sigma_{tot}$)''\n'r'$\sigma_{Qtl} = %.3f (\Delta = %.2f \sigma_{tot}$)'%(currentDistri.shape[0],100*(1-truncEvents*1./currentDistri.shape[0]),empiricMean,empiricStd,truncMean,(truncMean-empiricMean)/empiricStd, truncStd,(truncStd-empiricStd)/empiricStd, fitMean,(fitMean-empiricMean)/empiricStd,fitStd,(fitStd-empiricStd)/empiricStd,fwhmMean,(fwhmMean-empiricMean)/empiricStd,fwhmStd,(fwhmStd-empiricStd)/empiricStd,quantileMean,(quantileMean-empiricMean)/empiricStd,quantileStd,(quantileStd-empiricStd)/empiricStd),color = 'k',fontsize=16)

	plt.ylabel(r'$\mathrm{arbitrary\ units}$',fontsize = 22)


def weightedMean(data,weights):
	return np.average(data,weights=weights)

def weightedStd(data,weights):
	return np.sqrt(np.average((data-np.average(data,weights=weights))**2, weights=weights))

def weightedStdErr(data,weights):
	fourthMomentum = np.average((data-np.average(data,weights=weights))**4)
	n = data.shape[0]
	sig4 = np.average((data-np.average(data,weights=weights))**2, weights=weights)**2
	return (1./n*(fourthMomentum/sig4-3+2*n/(n-1))*sig4)



def get_weightedStdErr(data, dictPlot, method = 'Trunc'):
	if method == 'Empiric':
		#check if dataset contains datapoints
		if data.shape[0] == 0:
			return 0
		else:
			return weightedStdErr(data[:,0], data[:,1])
	else:
		centralData = data[((weightedMean(data[:,0],data[:,1])-4*weightedStd(data[:,0],data[:,1]))<data[:,0]) & (data[:,0] <(weightedMean(data[:,0],data[:,1])+4*weightedStd(data[:,0],data[:,1])))]
		if centralData.shape[0] == 0:
			return 0
		else:
			#(Empiric) mean on the truncated dataset
			if method == 'Trunc':
				return weightedStdErr(centralData[:,0], centralData[:,1])
			else:
				return 42


#calculates mean/std in the specified manner
def get_Quantity(data, dictPlot, quantity, method = 'Trunc'):
	#check if dataset contains datapoints
	if data.shape[0] == 0:
		return 0.01
	#Empiric mean over the whole data sample
	if method == 'Empiric':
		#Mean calculation
		if quantity == 'Mean':
			return weightedMean(data[:,0], data[:,1])
		#Standard deviation
		elif quantity == 'Std':
			return weightedStd(data[:,0], data[:,1])
		else:
			return 0.01
	else:
		#For the rest use truncated dataset, including all datapoints +/- 4 (empiric) sigma around the empiric mean 
		centralData = data[((weightedMean(data[:,0],data[:,1])-4*weightedStd(data[:,0],data[:,1]))<data[:,0]) & (data[:,0] <(weightedMean(data[:,0],data[:,1])+4*weightedStd(data[:,0],data[:,1])))]
		if method == 'Events':
			return centralData.shape[0]
		if centralData.shape[0] == 0:
			return 0.01
		else:
			#(Empiric) mean on the truncated dataset
			if method == 'Trunc':
				if quantity == 'Mean':
					return weightedMean(centralData[:,0], centralData[:,1])
				elif quantity == 'Std':
					return weightedStd(centralData[:,0], centralData[:,1])
				else:
					return 0.01
			#Fit Gaussian 
			else:
				num_bins = 50
				n, bins, patches = plt.hist(centralData[:,0], num_bins, facecolor='green', alpha=0.5, weights=centralData[:,1], normed=True)
				XCenter = (bins[:-1] + bins[1:])/2
				if method == 'Fit':
					if quantity == 'Mean':
						p0 = [1., 1., 1.]
						try:
							coeff, var_matrix = curve_fit(gauss,XCenter,n,p0=p0)
						except:
							coeff =[0,0.01,0.01]
						return coeff[1]
					elif quantity == 'Std':
						p0 = [1., 0., 15.]
						try:
							coeff, var_matrix = curve_fit(gauss,XCenter,n,p0=p0)
						except:
							coeff =[0,0.01,0.01]
						if not abs(coeff[2]) == 0:
							return abs(coeff[2])
						else:
							return 0.01
					else:
						return 0.01
				#Use Full width half maximum method with linear interpolation
				elif method == 'FWHM':
					FWHM_mean, FWHM_std = fwhm(XCenter,n)
					if quantity == 'Mean':
						return FWHM_mean
					elif quantity == 'Std':
						return FWHM_std/2.355
					else:
						return 0.01
				#Use Full width half maximum method with spline interpolation
				elif method == 'FWHMSpline':
					FWHM_mean, FWHM_std = fwhm(XCenter,n,'Spline')
					if quantity == 'Mean':
						return FWHM_mean
					elif quantity == 'Std':
						return FWHM_std/2.355
					else:
						return 0.01
				#Use quantile method
				elif method == 'Quantile':
					if quantity == 'Mean':
						return np.percentile(data[:,0], 50.)
					elif quantity == 'Std':
						return (np.percentile(data[:,0], 84.134)-np.percentile(data[:,0], 15.866))/2.
					else:
						return 0.01
				else:
					return 0.01

#Calculates next Point in which to calculate Alpha, based on 10% percent intervalls
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


#Perform alphaFit to extrapolate for perfect event topology
def alphaFit(dictPlot, alphaDataIn, bosonName='recoilslimmedMETs_Pt', mode='Response', saveName = 'dummy.pdf', method = 'AlphaInclInt'):
	XRange = np.zeros(9)
	YMean = np.zeros(9)
	YStd = np.zeros(9)

	#Data based on Trailing Jet Pt
	alphaData = alphaDataIn[(0<alphaDataIn[:,dictPlot['Jet1_Pt']])]
	try:
		minAlpha = getNextAlphaPoint(dictPlot, alphaData, 0)
		maxAlpha = getNextAlphaPoint(dictPlot, alphaData, minAlpha)
		for index in range(9):
			#Use all datapoints, including Trailing Jet Pt = 0 up to maxAlpha
			if method == 'AlphaInclInt':
				alphaDataLoop = alphaDataIn[(alphaDataIn[:,dictPlot['Jet1_Pt']]/alphaDataIn[:,dictPlot['Boson_Pt']]<maxAlpha)]
			#Use all datapoints, excluding Trailing Jet Pt = 0 up to maxAlpha
			elif method == 'AlphaExclInt':
				alphaDataLoop = alphaData[(alphaData[:,dictPlot['Jet1_Pt']]/alphaData[:,dictPlot['Boson_Pt']]<maxAlpha)]
			#Use only datapoints, from minAlpha to maxAlpha
			elif method == 'AlphaExcl':
				alphaDataLoop = alphaData[(alphaData[:,dictPlot['Jet1_Pt']]/alphaData[:,dictPlot['Boson_Pt']]>minAlpha) &(alphaData[:,dictPlot['Jet1_Pt']]/alphaData[:,dictPlot['Boson_Pt']]<maxAlpha)]
			#-U/ZPt
			if mode == 'Response':
				currentDistri = -alphaDataLoop[:,dictPlot[bosonName]]/alphaDataLoop[:,dictPlot['Boson_Pt']]
			#U+ZPt
			elif mode == 'Resolution':
				currentDistri = alphaDataLoop[:,dictPlot[bosonName]]+alphaDataLoop[:,dictPlot['Boson_Pt']]
			#fitDistri = currentDistri[((currentDistri.mean()-4*currentDistri.std())<currentDistri[:]) & (currentDistri[:] <(currentDistri.mean()+4*currentDistri.std()))]
			#XRange[index] = (maxAlpha+minAlpha)/2
			XRange[index] = maxAlpha
			#Calculate Mean and Std based on the truncated Mean method
			YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', method = 'Trunc')
			YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', method = 'Trunc')
			minAlpha = maxAlpha
			maxAlpha = getNextAlphaPoint(dictPlot, alphaData, minAlpha)
	except:
		p0 = [0.,1.]

	p0 = [0.,1.]

	plt.clf()
	#Perform linear Fit to datapoints to extrapolate for alpha = 0
	try:
		coeffMean, var_matrix = curve_fit(linear,XRange.transpose(),YMean.transpose(),p0=p0)
	except:
		coeffMean = [0.,0.]
	try:
		coeffStd, var_matrix = curve_fit(linear,XRange.transpose(),YStd.transpose(),p0=p0)
	except:
		coeffStd = [0.,0.]

	#Control Plots
	if mode == 'Response':
		plt.plot(XRange,YMean,'o')
		y = linear(XRange,*coeffMean)
		plt.plot(XRange,y)
		plt.ylabel(r'$<U_{\|\|} / p_T^H>$'' in GeV',fontsize = 20)
		plt.text(0.60*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.30*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$\mathrm{Events\,per\,Fit} = %i$''\n'r'$y = a \cdot x + b$''\n'r'$\mathrm{a} = %.2f$''\n'r'$\mathrm{b} = %.2f$''\n'%(alphaData.shape[0]/10,coeffMean[0], coeffMean[1]),color = 'k',fontsize=16)
	elif mode == 'Resolution':
		plt.plot(XRange,YStd,'o')
		y = linear(XRange,*coeffStd)
		plt.plot(XRange,y)
		plt.ylabel(r'$\sigma(<U_{\|\|} - p_T^H>)$'' in GeV',fontsize = 20)
		plt.text(0.60*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.30*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$\mathrm{Events\,per\,Fit} = %i$''\n'r'$y = a \cdot x + b$''\n'r'$\mathrm{a} = %.2f$''\n'r'$\mathrm{b} = %.2f$''\n'%(alphaData.shape[0]/10,coeffStd[0], coeffStd[1]),color = 'k',fontsize=16)
	plt.xlabel(r'$\alpha = p_{t}^{Jet1}/p_T^H}$',fontsize = 20)
	#if not saveName == 'dummy.pdf':
		#plt.savefig(saveName)
	plt.clf()

	#Return extrapolated Mean and Std
	return coeffMean[1], coeffStd[1]


#Full width half maximum method
def fwhm(x, y, method = 'Linear',k=10):
	"""
	Determine full-with-half-maximum of a peaked set of points, x and y.

	Assumes that there is only one peak present in the datasset.	The function
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

	#Take the points left and right of the maximum in case, more than 2 points were found
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
	#if too little points were found around the maximum
	elif len(roots) < 2:
		return 0, 0
	#if exactly two points were found
	else:
		return (roots[1] + roots[0])/2., abs(roots[1] - roots[0])


#gaussian function
def gauss(x, *p):
	A, mu, sigma = p
	return A*np.exp(-(x-mu)**2/(2.*sigma**2))

#linear function
def linear(x, *p):
	a, b = p
	return a*x + b

#Double gaussian function
def doublegauss(x, *p):
	A1, A2, mu1, mu2, sigma1, sigma2 = p
	return A1*np.exp(-(x-mu1)**2/(2.*sigma1**2)) + A2*np.exp(-(x-mu2)**2/(2.*sigma2**2))

#Add projection of MVA and PF Met on genMet
def add_MetProjections(config, inputDataPlot, dictPlot):
	#MVAMet
	if 'LongZCorrectedRecoil_MET' in dictPlot and 'LongZCorrectedRecoil_METPhi' in dictPlot and 'genMet_Pt' in dictPlot and 'genMet_Phi' in dictPlot:
		if not 'recoMetOnGenMetProjectionPar' in dictPlot:
			dictPlot['recoMetOnGenMetProjectionPar'] = inputDataPlot.shape[1]
			inputDataPlot = np.hstack((inputDataPlot, np.array(np.cos(inputDataPlot[:,dictPlot['genMet_Phi']]-inputDataPlot[:,dictPlot['LongZCorrectedRecoil_METPhi']])*(inputDataPlot[:,dictPlot['LongZCorrectedRecoil_MET']])).reshape(inputDataPlot.shape[0],1)))
		if not 'recoMetOnGenMetProjectionPerp' in dictPlot:
			dictPlot['recoMetOnGenMetProjectionPerp'] = inputDataPlot.shape[1]
			inputDataPlot = np.hstack((inputDataPlot, np.array(np.sin(inputDataPlot[:,dictPlot['genMet_Phi']]-inputDataPlot[:,dictPlot['LongZCorrectedRecoil_METPhi']])*(inputDataPlot[:,dictPlot['LongZCorrectedRecoil_MET']])).reshape(inputDataPlot.shape[0],1)))

	#Also for additional datasets
	for index in range(len(config['inputFile'])-1):
		if 'V%iLongZCorrectedRecoil_MET'%index in dictPlot and 'V%iLongZCorrectedRecoil_METPhi'%index in dictPlot and 'genMet_Pt' in dictPlot and 'genMet_Phi' in dictPlot:
			if not 'V%irecoMetOnGenMetProjectionPar'%index in dictPlot:
				dictPlot['V%irecoMetOnGenMetProjectionPar'%index] = inputDataPlot.shape[1]
				inputDataPlot = np.hstack((inputDataPlot, np.array(np.cos(inputDataPlot[:,dictPlot['genMet_Phi']]-inputDataPlot[:,dictPlot['V%iLongZCorrectedRecoil_METPhi'%index]])*(inputDataPlot[:,dictPlot['V%iLongZCorrectedRecoil_MET'%index]])).reshape(inputDataPlot.shape[0],1)))
			if not 'V%irecoMetOnGenMetProjectionPerp'%index in dictPlot:
				dictPlot['V%irecoMetOnGenMetProjectionPerp'%index] = inputDataPlot.shape[1]
				inputDataPlot = np.hstack((inputDataPlot, np.array(np.sin(inputDataPlot[:,dictPlot['genMet_Phi']]-inputDataPlot[:,dictPlot['V%iLongZCorrectedRecoil_METPhi'%index]])*(inputDataPlot[:,dictPlot['V%iLongZCorrectedRecoil_MET'%index]])).reshape(inputDataPlot.shape[0],1)))

	#Also for PF
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

	#Dictionary to transform ArtusVariables in Variables gained from MapAnalyzer
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


	#Read Input file
	filename = config['inputFile'][0]
	print('Loading dataset %s ...'%filename)

	if not os.path.exists(config['inputFile'][0].replace('.csv', '.pkl')): 
		print("Creating .pkl object for faster loading faster next time.")
		reader=csv.reader(open(filename,"rb"),delimiter=',') 
		datacsv=list(reader) 
		header = np.array(datacsv[0]).astype(np.str) 
		inputdatentot =np.array(datacsv[1:]).astype(np.float32) 
		pickle.dump([header, inputdatentot],open(config['inputFile'][0].replace('.csv', '.pkl'),'wb')) 
	else: 
		print('Found pickle format, using it instead of csv for faster loading process.')
		[header, inputdatentot] = pickle.load(open(config['inputFile'][0].replace('.csv', '.pkl'),'rb')) 
		print('Pickle loaded')



	"""
	reader=csv.reader(open(filename,"rb"),delimiter=',')
	datacsv=list(reader)
	header = np.array(datacsv[0]).astype(np.str)
	inputdatentot =np.array(datacsv[1:]).astype(np.float32)
	"""

	#Save input and replace name spaces from Artus
	dictInputTot = {}
	for index in range(0,header.shape[0]):
		if header[index] in ArtusDict:
			dictInputTot[ArtusDict[header[index]]] = index
		else:
			dictInputTot[header[index]] = index
				#print(header[index])

	#for name in dictInputTot:
			#print(name)

	#Variables to Plot
	plotnames = config[config['activePlotting']]['plotVariables']


	inputDataPlot = np.empty(shape=[inputdatentot.shape[0],0]).astype(np.float32)


	dictPlot = {}




	dt = int((time.time() - start))
	print('Elapsed time for loading dataset: ', dt)
	lastTime = time.time()

	#Add data to internal container
	for index, entry in enumerate(dictInputTot):
		if entry in plotnames:
			dictPlot[entry] = inputDataPlot.shape[1]
			inputDataPlot = np.hstack((inputDataPlot, np.array(inputdatentot[:,dictInputTot[entry]]).reshape(inputdatentot.shape[0],1)))

	#Load additional datasets
	if len(config['inputFile']) > 1:
		trainingheader =["LongZCorrectedRecoil_LongZ","LongZCorrectedRecoil_PerpZ","LongZCorrectedRecoil_Phi","PhiCorrectedRecoil_LongZ","PhiCorrectedRecoil_PerpZ","PhiCorrectedRecoil_Phi", "dmvamet_Pt", "dmvamet_Phi", "recoMetOnGenMetProjectionPar", "recoMetOnGenMetProjectionPerp","recoPfMetOnGenMetProjectionPar","recoPfMetOnGenMetProjectionPerp"]
		for index in range(len(config['inputFile'])-1):
			filename = config['inputFile'][index+1]
			print('Loading dataset %s ...'%filename)
			if not os.path.exists(config['inputFile'][0].replace('.csv', '.pkl')): 
				print("Creating .pkl object for faster loading faster next time.")
				reader=csv.reader(open(filename,"rb"),delimiter=',') 
				datacsv=list(reader) 
				header = np.array(datacsv[0]).astype(np.str) 
				inputdatentot =np.array(datacsv[1:]).astype(np.float32) 
				pickle.dump([header, inputdatentot],open(filename.replace('.csv', '.pkl'),'wb')) 
			else: 
				print('Found pickle format, using it instead of csv for faster loading process.')
				[header, inputdatentot] = pickle.load(open(filename.replace('.csv', '.pkl'),'rb')) 
				print('Pickle loaded')

			#take miminum dataset size for making stacking possible
			minSize = inputDataPlot.shape[0]
			if inputdatentot.shape[0] < minSize:
				minSize = inputdatentot.shape[0]
				print("Shrinked plotting data from %i to %i events to make simultaneous plotting possible"%(inputDataPlot.shape[0],minSize))
				inputDataPlot = inputDataPlot[:minSize,:]

			for indexHeader in range(0,header.shape[0]):
				if header[indexHeader] in trainingheader:
					dictPlot['V' + str(index)+header[indexHeader]] = inputDataPlot.shape[1]
					inputDataPlot = np.hstack((inputDataPlot, np.array(inputdatentot[:,indexHeader]).reshape(inputdatentot.shape[0],1)))
			dt = int(time.time() - lastTime)
			lastTime = time.time()
			print('Elapsed time for loading dataset: ', dt)


	#add projection on real Met for data from MapAnalyzer as it is not calculated there
	inputDataPlot, dictPlot = add_MetProjections(config, inputDataPlot, dictPlot)

	if 'trueBoson' in config:
		if config['trueBoson']:
			inputDataPlot[:,dictPlot['Boson_Pt']] = np.sqrt(inputDataPlot[:,dictPlot['Boson_Pt']]**2 + inputDataPlot[:,dictPlot['genMet_Pt']]**2 + 2*inputDataPlot[:,dictPlot['Boson_Pt']]*inputDataPlot[:,dictPlot['genMet_Pt']]*np.cos(inputDataPlot[:,dictPlot['genMet_Phi']]-inputDataPlot[:,dictPlot['Boson_Phi']]))
			inputDataPlot[:,dictPlot['Boson_Phi']] = inputDataPlot[:,dictPlot['Boson_Phi']] + np.arctan2(inputDataPlot[:,dictPlot['genMet_Pt']]*np.sin(inputDataPlot[:,dictPlot['genMet_Phi']]-inputDataPlot[:,dictPlot['Boson_Phi']]),inputDataPlot[:,dictPlot['Boson_Pt']]+inputDataPlot[:,dictPlot['genMet_Pt']]*np.cos(inputDataPlot[:,dictPlot['genMet_Phi']]-inputDataPlot[:,dictPlot['Boson_Phi']]))


	#if no event weights were specified
	if not 'eventWeight' in dictPlot:
		dictPlot['eventWeight'] = inputDataPlot.shape[1]
		inputDataPlot = np.hstack((inputDataPlot, np.ones((inputDataPlot.shape[0])).reshape(inputDataPlot.shape[0],1)))

	if 'LongZCorrectedRecoil_LongZ' in dictPlot:
		countFirst = inputDataPlot.shape[0]
		inputDataPlot = inputDataPlot[500>inputDataPlot[:,dictPlot['LongZCorrectedRecoil_LongZ']],:]
		if (countFirst - inputDataPlot.shape[0]) > 0:
			print(r'Deleted %i events exceeding a threshold of MVA recoil > 500 GeV'%(countFirst-inputDataPlot.shape[0]))

	print(dictPlot)
	#print(inputDataPlot.shape)



	dt = int((time.time() - start))
	print('Elapsed time for loading whole data: ', dt)

	return inputDataPlot, dictPlot



#Make shape plots for all variables used for plots
def make_Plot(variablename, inputData, dictPlot, outputdir):

	histData = inputData[:,dictPlot[variablename]]

	if histData.shape[0] == 0:
		return
	if not os.path.exists(outputdir):
		os.makedirs(outputdir)



	num_bins = 100

	if variablename == 'targetRecoilFromBDT' or variablename == 'targetRecoilFromSlimmed':
		n, bins, patches = plt.hist(histData, num_bins, facecolor='green', alpha=0.5, range=[-50, 50])
	else:
		n, bins, patches = plt.hist(histData, num_bins, facecolor='green', alpha=0.5)
	plt.xlabel(variablename,fontsize = 20)
	plt.ylabel('Hits',fontsize = 20)


	plt.tight_layout()
	plt.savefig((outputdir+variablename+".pdf"))
	plt.clf()
	return 0


def make_ResponseCorrectedPlot(config, XRange, YStd, YResponse, bosonName, targetvariable,	minrange,maxrange, stepwidth, ptmin,ptmax, labelname = 'MVA Met', relateVar = 'p_T^H', relateUnits = 'GeV'):

	plt.clf()
	ResCorr = YStd[:]/YResponse[:]
	binwidth = XRange[1:]-XRange[:-1]
	plt.plot(XRange[:-1]+stepwidth/2.,ResCorr,'o')
	plt.xlabel(targetvariable,fontsize = 20)
	plt.ylabel(r'Resolution / Response in GeV',fontsize = 20)
	if ptmax == 0:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/ResponseCorrected_%s_vs_%s.pdf" %(bosonName,targetvariable))
	else:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/ResponseCorrected%ito%iGeV_%s_vs_%s.pdf" %(ptmin,ptmax,bosonName,targetvariable))
	plt.figure(6)
	#plt.plot(XRange[:-1]+stepwidth/2.,ResCorr,'o-',label=labelname)
	plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),ResCorr[:],xerr=binwidth[:]/2.,fmt='o',label=labelname)
	#plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YStd[:],xerr=binwidth[:]/2.,yerr=YStdErr[:],fmt='o',label=labelname)
	#plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YStd[:],xerr=binwidth[:]/2.,fmt='o',label=labelname)
	plt.figure(0)
	plt.clf()



	return

def make_ResolutionPlot(config,plotData,dictPlot, bosonName, targetvariable,	minrange=42,maxrange=0, stepwidth=0, ptmin =0,ptmax=0, labelname = 'MVA Met', relateVar = 'p_T^H', relateUnits = 'GeV', binRanges = np.array([42])):


	if binRanges[0] == 42:
		if minrange == 42:
			minrange = plotData[:,dictPlot[targetvariable]].min()
		if maxrange == 0:
			maxrange = plotData[:,dictPlot[targetvariable]].max()
		if stepwidth == 0:
			stepwidth = (maxrange-minrange)/20.
		XRange = np.arange(minrange,maxrange,stepwidth)
		binwidth = np.zeros((XRange.shape[0]-1,1))
		binwidth[:] = stepwidth
	else:
		XRange = binRanges
		binwidth = binRanges[1:]-binRanges[:-1]

	YMean = np.zeros((XRange.shape[0]-1,1))
	YStd = np.zeros((XRange.shape[0]-1,1))
	YStdErr = np.zeros((XRange.shape[0]-1,1))

	print('Resolution %s versus %s'%(bosonName,targetvariable))
	#YValues
	for index in range(0,XRange.shape[0]-1):

		AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]
		currentDistri = (AlternativeDistri[:,dictPlot[bosonName]]+AlternativeDistri[:,dictPlot['Boson_Pt']]).reshape(AlternativeDistri.shape[0],1)
		currentDistri = np.hstack((currentDistri, np.array(AlternativeDistri[:,dictPlot['eventWeight']]).reshape(AlternativeDistri.shape[0],1)))
		#If no entries are assigned to this bin
		if currentDistri.shape[0] == 0:
			YMean[index] = 0
			YStd[index] = 0
		else:
			#if additional dataset is examined
			if bosonName[0] == 'V':
				#If several methods were chosen for the datasets
				if len(config['method']) > 1:
					if config['method'][int(bosonName[1])+1][0:5] == 'Alpha':
						YMean[index], YStd[index] = alphaFit(dictPlot, AlternativeDistri, bosonName,'Resolution', method = config['method'][int(bosonName[1])+1])
					else:
						YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][int(bosonName[1])+1])
						YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][int(bosonName[1])+1])
				#If one method is used for all datasets
				else:
					#If Alphamethod shall be used
					if config['method'][0][0:5] == 'Alpha':
						YMean[index], YStd[index] = alphaFit(dictPlot, AlternativeDistri, bosonName,'Resolution', method = config['method'][0])
					else:
						YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][0])
						YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][0])
			#First dataset
			else:
				if config['method'][0][0:5] == 'Alpha':
					YMean[index], YStd[index] = alphaFit(dictPlot, AlternativeDistri, bosonName,'Resolution', method = config['method'][0])
				else:
					YMean[index] = get_Quantity(currentDistri, dictPlot, 'Mean', config['method'][0])
					YStd[index] = get_Quantity(currentDistri, dictPlot, 'Std', config['method'][0])

			if config['studyVars']:
				quantityComparison(config, currentDistri, dictPlot, AlternativeDistri, bosonName)

			YStdErr[index] = get_weightedStdErr(currentDistri, dictPlot, config['method'][0])

			#plt.xlabel(r'$U_{||} - p_T^H$ at $%s = (%i - %i)\,\mathrm{%s}$'%(relateVar,XRange[index],XRange[index+1],relateUnits),fontsize = 20)
			plt.xlabel(r'$(U_{||} - p_T^H)\ \mathrm{in\ GeV}$',fontsize = 22)
			plt.text(0.04*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.9*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$(%i < %s < %i)\,\mathrm{%s}$'%(XRange[index],relateVar, XRange[index+1],relateUnits),color = 'k',fontsize=16, weight='bold', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
			plt.title(r'$\mathrm{Resolution}\ U_{||}\ - \ $ %s'%labelname[:-13]+'$', fontsize = 28, y = 1.02)
			if ptmax == 0:
				foldername = 'Resolution_%s_vs_%s' %(bosonName,targetvariable)
				if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
					os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
				bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
				os.system(bashCommand)
				plt.tight_layout()
				plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))
				"""
				if config['method'][0][0:5] == 'Alpha':
					alphaFit(dictPlot, AlternativeDistri, bosonName,'Resolution', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.pdf' %(foldername,foldername, index)), method = config['method'][0])
				else:
					alphaFit(dictPlot, AlternativeDistri, bosonName,'Resolution', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.pdf' %(foldername,foldername, index)))
				"""
			else:
				foldername = 'Resolution%ito%iGeV_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
				if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
					os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
				bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
				os.system(bashCommand)
				plt.tight_layout()
				plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))
				"""
				if config['method'][0][0:5] == 'Alpha':
					alphaFit(dictPlot, AlternativeDistri, bosonName,'Resolution', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.pdf' %(foldername,foldername, index)), method = config['method'][0])
				else:
					alphaFit(dictPlot, AlternativeDistri, bosonName,'Resolution', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.pdf' %(foldername,foldername, index)))
				"""
	plt.clf()
	plt.plot(XRange[:-1]+stepwidth/2.,YStd[:],'o-')
	plt.ylabel('(MET Boson PT_Long) - (True Boson Pt)',fontsize = 20)
	plt.xlabel(targetvariable,fontsize = 20)
	if ptmax == 0:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/Resolution_%s_vs_%s.pdf" %(bosonName,targetvariable))
	else:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/Resolution%ito%iGeV_%s_vs_%s.pdf" %(ptmin,ptmax,bosonName,targetvariable))
	plt.figure(5)
	#plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YStd[:],xerr=binwidth[:]/2.,yerr=YStdErr[:],fmt='o',label=labelname)
	plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YStd[:],xerr=binwidth[:]/2.,fmt='o',label=labelname)
	plt.figure(0)
	plt.clf()



	return XRange, YStd


def make_METResolutionPlot(config,plotData,dictPlot, bosonName, targetvariable,  minrange=42,maxrange=0, stepwidth=0, ptmin =0,ptmax=0, labelname = 'MVA Met', relateVar = 'p_T^H', relateUnits = 'GeV', binRanges = np.array([42])):


	#XRange = np.arange(plotData[:,targetindex].min(),plotData[:,targetindex].max(),(plotData[:,targetindex].max()-plotData[:,targetindex].min())/nbins)
	#if targetvariable == 'Boson_Pt':
		#binRanges=np.array([5,15,25,35,50,75,100,200])

	if binRanges[0] == 42:
		if minrange == 42:
			minrange = plotData[:,dictPlot[targetvariable]].min()
		if maxrange == 0:
			maxrange = plotData[:,dictPlot[targetvariable]].max()
		if stepwidth == 0:
			stepwidth = (maxrange-minrange)/20.
		XRange = np.arange(minrange,maxrange,stepwidth)
		binwidth = np.zeros((XRange.shape[0]-1,1))
		binwidth[:] = stepwidth
	else:
		XRange = binRanges
		binwidth = binRanges[1:]-binRanges[:-1]

	YMean = np.zeros((XRange.shape[0]-1,1))
	YStd = np.zeros((XRange.shape[0]-1,1))
	YStdErr = np.zeros((XRange.shape[0]-1,1))

	print('MET Resolution %s versus %s'%(bosonName,targetvariable))
	#YValues
	for index in range(0,XRange.shape[0]-1):

		AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]
		currentDistri = (AlternativeDistri[:,dictPlot[bosonName]]-AlternativeDistri[:,dictPlot['genMet_Pt']]).reshape(AlternativeDistri.shape[0],1)
		currentDistri = np.hstack((currentDistri, np.array(AlternativeDistri[:,dictPlot['eventWeight']]).reshape(AlternativeDistri.shape[0],1)))


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

		YStdErr[index] = get_weightedStdErr(currentDistri, dictPlot, config['method'][0])

		if config['studyVars']:
			quantityComparison(config, currentDistri, dictPlot, AlternativeDistri, bosonName)
		plt.text(0.04*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.9*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$(%i < %s < %i)\,\mathrm{%s}$'%(XRange[index],relateVar, XRange[index+1],relateUnits),color = 'k',fontsize=16, weight='bold', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
		plt.xlabel(r'$(E_T^{miss} - E_{T,true}^{miss})\ \mathrm{in\ GeV}$',fontsize = 22)
		plt.title(r'$\mathrm{Resolution}\ E_{T,||}^{miss}}\ - \ $  %s'%(labelname[:-13])+'$', fontsize = 28, y = 1.02)
		if ptmax == 0:
			foldername = 'METResolution_%s_vs_%s' %(bosonName,targetvariable)
			if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
				os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
				bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
				os.system(bashCommand)
			plt.tight_layout()
			plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))
		else:
			foldername = 'METResolution%ito%iGeV_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
			if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
				os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
				bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
				os.system(bashCommand)
			plt.tight_layout()
			plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))

	plt.clf()
	plt.plot(XRange[:-1]+stepwidth/2.,YStd[:],'o')
	plt.ylabel('(MET) - (gen MET)',fontsize = 20)
	plt.xlabel(targetvariable,fontsize = 20)
	if ptmax == 0:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/METResolution_%s_vs_%s.pdf" %(bosonName,targetvariable))
	else:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/METResolution%ito%iGeV_%s_vs_%s.pdf" %(ptmin,ptmax,bosonName,targetvariable))
	plt.figure(10)
	#print(YStdErr)
	#print(np.sqrt(YStdErr))
	#plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YStd[:],xerr=binwidth[:]/2.,yerr=YStdErr[:],fmt='o',label=labelname)
	plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YStd[:],xerr=binwidth[:]/2.,fmt='o',label=labelname)
	plt.figure(0)
	plt.clf()


	return 0

def make_METResolutionPerpPlot(config,plotData,dictPlot, bosonName, targetvariable,  minrange=42,maxrange=0, stepwidth=0, ptmin =0,ptmax=0, labelname = 'MVA Met', relateVar = 'p_T^H', relateUnits = 'GeV', binRanges = np.array([42])):


	#XRange = np.arange(plotData[:,targetindex].min(),plotData[:,targetindex].max(),(plotData[:,targetindex].max()-plotData[:,targetindex].min())/nbins)
	#if targetvariable == 'Boson_Pt':
		#binRanges=np.array([5,15,25,35,50,75,100,200])

	if binRanges[0] == 42:
		if minrange == 42:
			minrange = plotData[:,dictPlot[targetvariable]].min()
		if maxrange == 0:
			maxrange = plotData[:,dictPlot[targetvariable]].max()
		if stepwidth == 0:
			stepwidth = (maxrange-minrange)/20.
		XRange = np.arange(minrange,maxrange,stepwidth)
		binwidth = np.zeros((XRange.shape[0]-1,1))
		binwidth[:] = stepwidth
	else:
		XRange = binRanges
		binwidth = binRanges[1:]-binRanges[:-1]

	YMean = np.zeros((XRange.shape[0]-1,1))
	YStd = np.zeros((XRange.shape[0]-1,1))
	YStdErr = np.zeros((XRange.shape[0]-1,1))

	print('MET Resolution Perp %s versus %s'%(bosonName,targetvariable))
	#YValues
	for index in range(0,XRange.shape[0]-1):

		AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]
		currentDistri = (AlternativeDistri[:,dictPlot[bosonName]]).reshape(AlternativeDistri.shape[0],1)
		currentDistri = np.hstack((currentDistri, np.array(AlternativeDistri[:,dictPlot['eventWeight']]).reshape(AlternativeDistri.shape[0],1)))


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

		YStdErr[index] = get_weightedStdErr(currentDistri, dictPlot, config['method'][0])

		if config['studyVars']:
			quantityComparison(config, currentDistri, dictPlot, AlternativeDistri, bosonName)
		plt.xlabel(r'$E_{T,\bot}^{miss}\ \mathrm{in\ GeV}$',fontsize = 22)
		plt.text(0.04*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.9*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$(%i < %s < %i)\,\mathrm{%s}$'%(XRange[index],relateVar, XRange[index+1],relateUnits),color = 'k',fontsize=16, weight='bold', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
		plt.title(r'$\mathrm{Resolution}\ E_{T,bot}^{miss}\ - \ $ %s'%(labelname[:-13])+'$', fontsize = 28, y = 1.02)
		if ptmax == 0:
			foldername = 'METResolutionPerp_%s_vs_%s' %(bosonName,targetvariable)
			if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
				os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
				bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
				os.system(bashCommand)
			plt.tight_layout()
			plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))
		else:
			foldername = 'METResolutionPerp%ito%iGeV_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
			if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
				os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
				bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
				os.system(bashCommand)
			plt.tight_layout()
			plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))

	plt.clf()
	plt.plot(XRange[:-1]+stepwidth/2.,YStd[:],'o')
	plt.ylabel('(MET) - (gen MET)',fontsize = 20)
	plt.xlabel(targetvariable,fontsize = 20)
	if ptmax == 0:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/METResolutionPerp_%s_vs_%s.pdf" %(bosonName,targetvariable))
	else:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/METResolutionPerp%ito%iGeV_%s_vs_%s.pdf" %(ptmin,ptmax,bosonName,targetvariable))
	plt.figure(11)
	#print(YStdErr)
	#plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YStd[:],xerr=binwidth[:]/2.,yerr=YStdErr[:],fmt='o',label=labelname)
	plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YStd[:],xerr=binwidth[:]/2.,fmt='o',label=labelname)
	plt.figure(0)
	plt.clf()


	return 0



def make_ResponsePlot(config, plotData,dictPlot, bosonName, targetvariable,  minrange=42,maxrange=0, stepwidth=0, ptmin=0, ptmax=0, labelname = 'MVA Met', relateVar = 'p_T^H', relateUnits = 'GeV', binRanges = np.array([42])):

	if binRanges[0] == 42:
		if minrange == 42:
			minrange = plotData[:,dictPlot[targetvariable]].min()
		if maxrange == 0:
			maxrange = plotData[:,dictPlot[targetvariable]].max()
		if stepwidth == 0:
			stepwidth = (maxrange-minrange)/20.
		XRange = np.arange(minrange,maxrange,stepwidth)
		binwidth = np.zeros((XRange.shape[0]-1,1))
		binwidth[:] = stepwidth
	else:
		XRange = binRanges
		binwidth = binRanges[1:]-binRanges[:-1]

	YMean = np.zeros((XRange.shape[0]-1,1))
	YStd = np.zeros((XRange.shape[0]-1,1))
	YStdErr = np.zeros((XRange.shape[0]-1,1))
	print('Response %s versus %s'%(bosonName,targetvariable))



	#YValues
	ignoredEntries = 0
	for index in range(0,XRange.shape[0]-1):

		AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]

		#Inverse definition
		#currentDistri = (-AlternativeDistri[:,dictPlot['Boson_Pt']]/AlternativeDistri[:,dictPlot[bosonName]]).reshape(AlternativeDistri.shape[0],1)
		currentDistri = (-AlternativeDistri[:,dictPlot[bosonName]]/AlternativeDistri[:,dictPlot['Boson_Pt']]).reshape(AlternativeDistri.shape[0],1)
		currentDistri = np.hstack((currentDistri, np.array(AlternativeDistri[:,dictPlot['eventWeight']]).reshape(AlternativeDistri.shape[0],1)))
		if currentDistri.shape[0] == 0:
			YMean[index] = 0
			YStd[index] = 0
		else:

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

			YStdErr[index] = get_weightedStdErr(currentDistri, dictPlot, config['method'][0])
			if config['studyVars']:
				quantityComparison(config, currentDistri, dictPlot, AlternativeDistri, bosonName)
			plt.xlabel(r'$U_{||} / -p_T^H$',fontsize = 22)
			plt.text(0.04*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.9*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$(%i < %s < %i)\,\mathrm{%s}$'%(XRange[index],relateVar, XRange[index+1],relateUnits),color = 'k',fontsize=16, weight='bold', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
			plt.title(r'$\mathrm{Response}\ U_{||}-$ %s'%labelname[:-13]+'$', fontsize = 28, y = 1.02)
			if ptmax == 0:
				foldername = 'Response_%s_vs_%s' %(bosonName,targetvariable)
				if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
					os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
					bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
					os.system(bashCommand)
				plt.tight_layout()
				plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))
				"""
				if config['method'][0][0:5] == 'Alpha':
					alphaFit(dictPlot, AlternativeDistri, bosonName,'Response', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.pdf' %(foldername,foldername, index)), method = config['method'][0])
				else:
					alphaFit(dictPlot, AlternativeDistri, bosonName,'Response', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.pdf' %(foldername,foldername, index)))
				"""
			else:
				foldername = 'Response%ito%iGeV_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
				if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
					os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
					bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
					os.system(bashCommand)
				plt.tight_layout()
				plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))

				"""
				if config['method'][0][0:5] == 'Alpha':
					alphaFit(dictPlot, AlternativeDistri, bosonName,'Response', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.pdf' %(foldername,foldername, index)), method = config['method'][0])
				else:
					alphaFit(dictPlot, AlternativeDistri, bosonName,'Response', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.pdf' %(foldername,foldername, index)))
				"""

	plt.clf()
	plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o-')

	plt.xlabel(targetvariable,fontsize = 20)
	if ptmax == 0:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/Response_%s_vs_%s.pdf" %(bosonName,targetvariable))
	else:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/Response%ito%iGeV_%s_vs_%s.pdf" %(ptmin,ptmax,bosonName,targetvariable))
	plt.clf()
	plt.figure(4)
	#plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YStd[:],xerr=binwidth[:]/2.,yerr=YStdErr[:],fmt='o',label=labelname)
	plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YMean[:],xerr=binwidth[:]/2.,fmt='o',label=labelname)
	#print(YMean)
	#print(YStd)
	plt.figure(0)


	return YMean

def make_ResponseInvertedPlot(config, plotData,dictPlot, bosonName, targetvariable,  minrange=42,maxrange=0, stepwidth=0, ptmin=0, ptmax=0, labelname = 'MVA Met', relateVar = 'p_T^H', relateUnits = 'GeV', binRanges = np.array([42])):

	if binRanges[0] == 42:
		if minrange == 42:
			minrange = plotData[:,dictPlot[targetvariable]].min()
		if maxrange == 0:
			maxrange = plotData[:,dictPlot[targetvariable]].max()
		if stepwidth == 0:
			stepwidth = (maxrange-minrange)/20.
		XRange = np.arange(minrange,maxrange,stepwidth)
		binwidth = np.zeros((XRange.shape[0]-1,1))
		binwidth[:] = stepwidth
	else:
		XRange = binRanges
		binwidth = binRanges[1:]-binRanges[:-1]

	YMean = np.zeros((XRange.shape[0]-1,1))
	YStd = np.zeros((XRange.shape[0]-1,1))
	YStdErr = np.zeros((XRange.shape[0]-1,1))
	print('Response Inverted %s versus %s'%(bosonName,targetvariable))



	#YValues
	ignoredEntries = 0
	for index in range(0,XRange.shape[0]-1):

		AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]

		#Inverse definition
		currentDistri = (-AlternativeDistri[:,dictPlot['Boson_Pt']]/AlternativeDistri[:,dictPlot[bosonName]]).reshape(AlternativeDistri.shape[0],1)
		#currentDistri = (-AlternativeDistri[:,dictPlot[bosonName]]/AlternativeDistri[:,dictPlot['Boson_Pt']]).reshape(AlternativeDistri.shape[0],1)
		currentDistri = np.hstack((currentDistri, np.array(AlternativeDistri[:,dictPlot['eventWeight']]).reshape(AlternativeDistri.shape[0],1)))
		if currentDistri.shape[0] == 0:
			YMean[index] = 0
			YStd[index] = 0
		else:

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

			YStdErr[index] = get_weightedStdErr(currentDistri, dictPlot, config['method'][0])
			if config['studyVars']:
				quantityComparison(config, currentDistri, dictPlot, AlternativeDistri, bosonName)
			plt.xlabel(r'$ -p_T^H/ U_{||}$',fontsize = 22)
			plt.text(0.04*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.9*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$(%i < %s < %i)\,\mathrm{%s}$'%(XRange[index],relateVar, XRange[index+1],relateUnits),color = 'k',fontsize=16, weight='bold', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
			plt.title(r'$\mathrm{Response^{-1}}\ U_{||}\ -\ $ %s'%labelname[:-13]+'$', fontsize = 28, y = 1.02)
			if ptmax == 0:
				foldername = 'ResponseInverted_%s_vs_%s' %(bosonName,targetvariable)
				if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
					os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
					bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
					os.system(bashCommand)
				plt.tight_layout()
				plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))
				"""
				if config['method'][0][0:5] == 'Alpha':
					alphaFit(dictPlot, AlternativeDistri, bosonName,'Response', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.pdf' %(foldername,foldername, index)), method = config['method'][0])
				else:
					alphaFit(dictPlot, AlternativeDistri, bosonName,'Response', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.pdf' %(foldername,foldername, index)))
				"""
			else:
				foldername = 'ResponseInverted%ito%iGeV_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
				if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
					os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
					bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
					os.system(bashCommand)
				plt.tight_layout()
				plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))

				"""
				if config['method'][0][0:5] == 'Alpha':
					alphaFit(dictPlot, AlternativeDistri, bosonName,'Response', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.pdf' %(foldername,foldername, index)), method = config['method'][0])
				else:
					alphaFit(dictPlot, AlternativeDistri, bosonName,'Response', (config['outputDir'] + 'ControlPlots/SingleDistributions/%s/alpha_%s_%i.pdf' %(foldername,foldername, index)))
				"""

	plt.clf()
	plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o-')

	plt.xlabel(targetvariable,fontsize = 20)
	if ptmax == 0:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/ResponseInv_%s_vs_%s.pdf" %(bosonName,targetvariable))
	else:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/ResponseInv%ito%iGeV_%s_vs_%s.pdf" %(ptmin,ptmax,bosonName,targetvariable))
	plt.clf()
	plt.figure(12)
	#plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YStd[:],xerr=binwidth[:]/2.,yerr=YStdErr[:],fmt='o',label=labelname)
	plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YMean[:],xerr=binwidth[:]/2.,fmt='o',label=labelname)
	#print(YMean)
	#print(YStd)
	plt.figure(0)

	return 0

def make_METResponseInvertedPlot(config, plotData,dictPlot, bosonName, targetvariable,	minrange=42,maxrange=0, stepwidth=0, ptmin=0, ptmax=0, labelname = 'MVA Met', relateVar = 'p_T^H', relateUnits = 'GeV', binRanges = np.array([42])):

	if binRanges[0] == 42:
		if minrange == 42:
			minrange = plotData[:,dictPlot[targetvariable]].min()
		if maxrange == 0:
			maxrange = plotData[:,dictPlot[targetvariable]].max()
		if stepwidth == 0:
			stepwidth = (maxrange-minrange)/20.
		XRange = np.arange(minrange,maxrange,stepwidth)
		binwidth = np.zeros((XRange.shape[0]-1,1))
		binwidth[:] = stepwidth
	else:
		XRange = binRanges
		binwidth = binRanges[1:]-binRanges[:-1]

	YMean = np.zeros((XRange.shape[0]-1,1))
	YStd = np.zeros((XRange.shape[0]-1,1))
	YStdErr = np.zeros((XRange.shape[0]-1,1))

	print('MET Response Inverted %s versus %s'%(bosonName,targetvariable))

	#YValues
	ignoredEntries = 0
	for index in range(0,XRange.shape[0]-1):

		AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]

		currentDistri = (AlternativeDistri[:,dictPlot['genMet_Pt']]/AlternativeDistri[:,dictPlot[bosonName]]).reshape(AlternativeDistri.shape[0],1)
		currentDistri = np.hstack((currentDistri, np.array(AlternativeDistri[:,dictPlot['eventWeight']]).reshape(AlternativeDistri.shape[0],1)))

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

		YStdErr[index] = get_weightedStdErr(currentDistri, dictPlot, config['method'][0])
		if config['studyVars']:
			quantityComparison(config, currentDistri, dictPlot, AlternativeDistri, bosonName)
		plt.xlabel(r'$E_{T,gen}^{miss}/E_{T}^{miss}$',fontsize = 22)
		plt.text(0.04*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.9*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$(%i < %s < %i)\,\mathrm{%s}$'%(XRange[index],relateVar, XRange[index+1],relateUnits),color = 'k',fontsize=16, weight='bold', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
		plt.title(r'$\mathrm{Response^{-1}}\ E_{T,||}{miss}\ -\ $ %s'%labelname[:-13]+'$', fontsize = 28, y = 1.02)
		if ptmax == 0:
			foldername = 'METResponseInverted_%s_vs_%s' %(bosonName,targetvariable)
			if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
				os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
				bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
				os.system(bashCommand)
			plt.tight_layout()
			plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))
		else:
			foldername = 'METResponseInverted%ito%iGeV_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
			if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
				os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
				bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
				os.system(bashCommand)
			plt.tight_layout()
			plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))

	plt.clf()
	plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o')

	plt.xlabel(targetvariable,fontsize = 20)
	if ptmax == 0:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/METResponseInverted_%s_vs_%s.pdf" %(bosonName,targetvariable))
	else:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/METResponseInverted%ito%iGeV_%s_vs_%s.pdf" %(ptmin,ptmax,bosonName,targetvariable))
	plt.clf()
	plt.figure(13)
	#plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YStd[:],xerr=binwidth[:]/2.,yerr=YStdErr[:],fmt='o',label=labelname)
	plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YMean[:],xerr=binwidth[:]/2.,fmt='o',label=labelname)

	plt.figure(0)


	return 0

def make_METResponsePlot(config, plotData,dictPlot, bosonName, targetvariable,	minrange=42,maxrange=0, stepwidth=0, ptmin=0, ptmax=0, labelname = 'MVA Met', relateVar = 'p_T^H', relateUnits = 'GeV', binRanges = np.array([42])):

	if binRanges[0] == 42:
		if minrange == 42:
			minrange = plotData[:,dictPlot[targetvariable]].min()
		if maxrange == 0:
			maxrange = plotData[:,dictPlot[targetvariable]].max()
		if stepwidth == 0:
			stepwidth = (maxrange-minrange)/20.
		XRange = np.arange(minrange,maxrange,stepwidth)
		binwidth = np.zeros((XRange.shape[0]-1,1))
		binwidth[:] = stepwidth
	else:
		XRange = binRanges
		binwidth = binRanges[1:]-binRanges[:-1]

	YMean = np.zeros((XRange.shape[0]-1,1))
	YStd = np.zeros((XRange.shape[0]-1,1))
	YStdErr = np.zeros((XRange.shape[0]-1,1))

	print('MET Response %s versus %s'%(bosonName,targetvariable))

	#YValues
	ignoredEntries = 0
	for index in range(0,XRange.shape[0]-1):

		AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]

		currentDistri = (AlternativeDistri[:,dictPlot[bosonName]]/AlternativeDistri[:,dictPlot['genMet_Pt']]).reshape(AlternativeDistri.shape[0],1)
		currentDistri = np.hstack((currentDistri, np.array(AlternativeDistri[:,dictPlot['eventWeight']]).reshape(AlternativeDistri.shape[0],1)))

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

		YStdErr[index] = get_weightedStdErr(currentDistri, dictPlot, config['method'][0])
		if config['studyVars']:
			quantityComparison(config, currentDistri, dictPlot, AlternativeDistri, bosonName)
		plt.xlabel(r'$E_{t}^{miss}  / E_{t,gen}^{miss}$', fontsize = 22)
		plt.text(0.04*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.9*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$(%i < %s < %i)\,\mathrm{%s}$'%(XRange[index],relateVar, XRange[index+1],relateUnits),color = 'k',fontsize=16, weight='bold', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
		plt.title(r'$\mathrm{Response}\ E_{T,||}^{miss}\ -\ $ %s'%labelname[:-13]+'$', fontsize = 28, y = 1.02)
		if ptmax == 0:
			foldername = 'METResponse_%s_vs_%s' %(bosonName,targetvariable)
			if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
				os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
				bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
				os.system(bashCommand)
			plt.tight_layout()
			plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))
		else:
			foldername = 'METResponse%ito%iGeV_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
			if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
				os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
				bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
				os.system(bashCommand)
			plt.tight_layout()
			plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))

	plt.clf()
	plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o')

	plt.xlabel(targetvariable,fontsize = 20)
	if ptmax == 0:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/METResponse_%s_vs_%s.pdf" %(bosonName,targetvariable))
	else:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/METResponse%ito%iGeV_%s_vs_%s.pdf" %(ptmin,ptmax,bosonName,targetvariable))
	plt.clf()
	plt.figure(9)
	#plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YStd[:],xerr=binwidth[:]/2.,yerr=YStdErr[:],fmt='o',label=labelname)
	plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YMean[:],xerr=binwidth[:]/2.,fmt='o',label=labelname)

	plt.figure(0)


	return 0



def make_ResolutionPerpPlot(config,plotData,dictPlot, bosonName, targetvariable,	minrange=42,maxrange=0, stepwidth=0, ptmin =0,ptmax=0, labelname = 'MVA Met', relateVar = 'p_T^H', relateUnits = 'GeV', binRanges = np.array([42])):

	if binRanges[0] == 42:
		if minrange == 42:
			minrange = plotData[:,dictPlot[targetvariable]].min()
		if maxrange == 0:
			maxrange = plotData[:,dictPlot[targetvariable]].max()
		if stepwidth == 0:
			stepwidth = (maxrange-minrange)/20.
		XRange = np.arange(minrange,maxrange,stepwidth)
		binwidth = np.zeros((XRange.shape[0]-1,1))
		binwidth[:] = stepwidth
	else:
		XRange = binRanges
		binwidth = binRanges[1:]-binRanges[:-1]

	YMean = np.zeros((XRange.shape[0]-1,1))
	YStd = np.zeros((XRange.shape[0]-1,1))
	YStdErr = np.zeros((XRange.shape[0]-1,1))


	print('Resolution Perp %s versus %s'%(bosonName,targetvariable))
	#YValues
	for index in range(0,XRange.shape[0]-1):

		AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]
		currentDistri = AlternativeDistri[:,dictPlot[bosonName]].reshape(AlternativeDistri.shape[0],1)
		currentDistri = np.hstack((currentDistri, np.array(AlternativeDistri[:,dictPlot['eventWeight']]).reshape(AlternativeDistri.shape[0],1)))


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

		YStdErr[index] = get_weightedStdErr(currentDistri, dictPlot, config['method'][0])
		if config['studyVars']:
			quantityComparison(config, currentDistri, dictPlot, AlternativeDistri, bosonName)
		plt.xlabel(r'$U_\bot\ \mathrm{in\ GeV}$',fontsize = 22)
		plt.text(0.04*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.9*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$(%i < %s < %i)\,\mathrm{%s}$'%(XRange[index],relateVar, XRange[index+1],relateUnits),color = 'k',fontsize=16, weight='bold', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
		plt.title(r'$\mathrm{Resolution}\ U_{\bot}\ -\ $ %s'%labelname[:-13]+'$', fontsize = 28, y = 1.02)
		if index < 12:
			if ptmax == 0:
				foldername = 'ResolutionPerp_%s_vs_%s' %(bosonName,targetvariable)
				if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
					os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
					bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
					os.system(bashCommand)
				plt.tight_layout()
				plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))
			else:
				foldername = 'ResolutionPerp%ito%iGeV_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
				if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
					os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
					bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
					os.system(bashCommand)
				plt.tight_layout()
				plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))

	plt.clf()
	plt.plot(XRange[:-1]+stepwidth/2.,YStd[:],'o')
	plt.ylabel('MET Boson PT_Perp',fontsize = 20)
	plt.xlabel(targetvariable,fontsize = 20)
	if ptmax == 0:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/ResolutionPerp_%s_vs_%s.pdf" %(bosonName,targetvariable))
	else:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/ResolutionPerp%ito%iGeV_%s_vs_%s.pdf" %(ptmin,ptmax,bosonName,targetvariable))

	plt.figure(8)
	#plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YStd[:],xerr=binwidth[:]/2.,yerr=YStdErr[:],fmt='o',label=labelname)
	plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YStd[:],xerr=binwidth[:]/2.,fmt='o',label=labelname)
	plt.figure(0)
	plt.clf()


	return

def make_ResponsePerpPlot(config, plotData,dictPlot, bosonName, targetvariable,  minrange=42,maxrange=0, stepwidth=0, ptmin=0, ptmax=0, labelname = 'MVA Met', relateVar = 'p_T^H', relateUnits = 'GeV', binRanges = np.array([42])):

	if binRanges[0] == 42:
		if minrange == 42:
			minrange = plotData[:,dictPlot[targetvariable]].min()
		if maxrange == 0:
			maxrange = plotData[:,dictPlot[targetvariable]].max()
		if stepwidth == 0:
			stepwidth = (maxrange-minrange)/20.
		XRange = np.arange(minrange,maxrange,stepwidth)
		binwidth = np.zeros((XRange.shape[0]-1,1))
		binwidth[:] = stepwidth
	else:
		XRange = binRanges
		binwidth = binRanges[1:]-binRanges[:-1]

	YMean = np.zeros((XRange.shape[0]-1,1))
	YStd = np.zeros((XRange.shape[0]-1,1))
	YStdErr = np.zeros((XRange.shape[0]-1,1))

	print('Response Perp %s versus %s'%(bosonName,targetvariable))


	#YValues
	ignoredEntries = 0
	for index in range(0,XRange.shape[0]-1):

		AlternativeDistri = plotData[(XRange[index]<plotData[:,dictPlot[targetvariable]]) & (XRange[index+1]>plotData[:,dictPlot[targetvariable]])]

		currentDistri = AlternativeDistri[:,dictPlot[bosonName]].reshape(AlternativeDistri.shape[0],1)
		currentDistri = np.hstack((currentDistri, np.array(AlternativeDistri[:,dictPlot['eventWeight']]).reshape(AlternativeDistri.shape[0],1)))

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

		YStdErr[index] = get_weightedStdErr(currentDistri, dictPlot, config['method'][0])
		if config['studyVars']:
			quantityComparison(config, currentDistri, dictPlot, AlternativeDistri, bosonName)
		plt.xlabel(r'$U_\bot$ in GeV',fontsize = 20)
		plt.text(0.04*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.9*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$(%i < %s < %i)\,\mathrm{%s}$'%(XRange[index],relateVar, XRange[index+1],relateUnits),color = 'k',fontsize=16, weight='bold', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
		plt.title('Response Perp %s'%labelname[:-13]+'$', fontsize = 20)
		if ptmax == 0:
			foldername = 'ResponsePerp_%s_vs_%s' %(bosonName,targetvariable)
			if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
				os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
				bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
				os.system(bashCommand)
			plt.tight_layout()
			plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))
		else:
			foldername = 'ResponsePerp%ito%iGeV_%s_vs_%s' %(ptmin,ptmax,bosonName,targetvariable)
			if not os.path.exists(config['outputDir'] + 'ControlPlots/SingleDistributions/%s/'%foldername):
				os.makedirs(config['outputDir'] + 'ControlPlots/SingleDistributions/%s'%foldername)
				bashCommand = 'cp index.php %sControlPlots/SingleDistributions/%s/'%(config['outputDir'],foldername)
				os.system(bashCommand)
			plt.tight_layout()
			plt.savefig((config['outputDir'] + 'ControlPlots/SingleDistributions/%s/%s_%i.pdf' %(foldername,foldername, index)))
	plt.clf()
	plt.plot(XRange[:-1]+stepwidth/2.,YMean[:],'o')

	plt.xlabel(targetvariable,fontsize = 20)
	if ptmax == 0:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/ResponsePerp_%s_vs_%s.pdf" %(bosonName,targetvariable))
	else:
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "ControlPlots/ResponsePerp%ito%iGeV_%s_vs_%s.pdf" %(ptmin,ptmax,bosonName,targetvariable))
	plt.clf()
	plt.figure(7)
	#plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YStd[:],xerr=binwidth[:]/2.,yerr=YStdErr[:],fmt='o',label=labelname)
	plt.errorbar((XRange[:-1]+binwidth[:].transpose()/2.).transpose(),YMean[:],xerr=binwidth[:]/2.,fmt='o',label=labelname)
	plt.figure(0)



	return



def make_ControlPlots(config, plotData,dictPlot, bosonName, targetvariable,  minrange=42,maxrange=0, stepwidth=0, ptmin=0,ptmax=0, labelname = 'MVA Met', relateVar = 'p_T^H', relateUnits = 'GeV', binRanges = np.array([42])):

	bosonNameLong = bosonName + '_LongZ'
	bosonNamePerp = bosonName + '_PerpZ'
	maxrange += stepwidth
	if not os.path.exists((config['outputDir'] + 'ControlPlots/SingleDistributions/')):
		os.makedirs((config['outputDir'] + 'ControlPlots/SingleDistributions/'))

	if bosonNameLong in dictPlot:
		XRange, YVariance = make_ResolutionPlot(config, plotData, dictPlot, bosonNameLong, targetvariable,	minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits, binRanges = binRanges)
		YResponse = make_ResponsePlot(config, plotData, dictPlot, bosonNameLong, targetvariable,	minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits, binRanges = binRanges)
		make_ResponseInvertedPlot(config, plotData, dictPlot, bosonNameLong, targetvariable,	minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits, binRanges = binRanges)
		make_ResponseCorrectedPlot(config, XRange, YVariance, YResponse, bosonNameLong, targetvariable,  minrange,maxrange, stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits)
	if bosonNamePerp in dictPlot:
		make_ResolutionPerpPlot(config, plotData, dictPlot, bosonNamePerp, targetvariable,	minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits, binRanges = binRanges)
		make_ResponsePerpPlot(config, plotData, dictPlot, bosonNamePerp, targetvariable,	minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits, binRanges = binRanges)

	if bosonName == "LongZCorrectedRecoil" and "recoMetOnGenMetProjectionPar" in dictPlot and not plotData[:,dictPlot["recoMetOnGenMetProjectionPar"]].mean() == -999:
		make_METResponsePlot(config, plotData, dictPlot, "recoMetOnGenMetProjectionPar", targetvariable,	minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits, binRanges = binRanges)
		make_METResponseInvertedPlot(config, plotData, dictPlot, "recoMetOnGenMetProjectionPar", targetvariable,	minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits, binRanges = binRanges)
		make_METResolutionPlot(config, plotData, dictPlot, "recoMetOnGenMetProjectionPar", targetvariable,	minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits, binRanges = binRanges)
		make_METResolutionPerpPlot(config, plotData, dictPlot, "recoMetOnGenMetProjectionPerp", targetvariable,	minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits, binRanges = binRanges)

	if bosonName[0] == "V" and "V%srecoMetOnGenMetProjectionPar"%bosonName[1] in dictPlot and not plotData[:,dictPlot["V%srecoMetOnGenMetProjectionPar"%bosonName[1]]].mean() == -999:
		make_METResponsePlot(config, plotData, dictPlot, "V%srecoMetOnGenMetProjectionPar"%bosonName[1], targetvariable,	minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits, binRanges = binRanges)
		make_METResponseInvertedPlot(config, plotData, dictPlot, "V%srecoMetOnGenMetProjectionPar"%bosonName[1], targetvariable,	minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits, binRanges = binRanges)
		make_METResolutionPlot(config, plotData, dictPlot, "V%srecoMetOnGenMetProjectionPar"%bosonName[1], targetvariable,	minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits, binRanges = binRanges)
		make_METResolutionPerpPlot(config, plotData, dictPlot, "V%srecoMetOnGenMetProjectionPerp"%bosonName[1], targetvariable,	minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits, binRanges = binRanges)

	if bosonName == "recoilslimmedMETs" and "recoPfMetOnGenMetProjectionPar" in dictPlot and not plotData[:,dictPlot["recoPfMetOnGenMetProjectionPar"]].mean() == -999:
		make_METResponsePlot(config, plotData, dictPlot, "recoPfMetOnGenMetProjectionPar", targetvariable,	minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits, binRanges = binRanges)
		make_METResponseInvertedPlot(config, plotData, dictPlot, "recoPfMetOnGenMetProjectionPar", targetvariable,	minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits, binRanges = binRanges)
		make_METResolutionPlot(config, plotData, dictPlot, "recoPfMetOnGenMetProjectionPar", targetvariable,	minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits, binRanges = binRanges)
		make_METResolutionPerpPlot(config, plotData, dictPlot, "recoPfMetOnGenMetProjectionPerp", targetvariable,	minrange,maxrange,stepwidth, ptmin, ptmax, labelname, relateVar, relateUnits, binRanges = binRanges)



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

		YResponse = make_ResponsePlot(config, lowerData, dictPlot, bosonNameLong, 'Boson_Pt',  10,200,10,0,0,'PF Met','p_T^H','GeV')
		XRange, YVariance = make_ResolutionPlot(config, lowerData, dictPlot, bosonNameLong, 'Boson_Pt',  10,200,10,0,0,'PF Met','p_T^H','GeV')

	if 'recoilslimmedMETs_LongZ' in dictPlot:
		bosonNameLong = 'recoilslimmedMETs_LongZ'
		bosonNamePerp = 'recoilslimmedMETs_PerpZ'

		if not os.path.exists((config['outputDir'] + 'ControlPlots/SingleDistributions/')):
			os.makedirs((config['outputDir'] + 'ControlPlots/SingleDistributions/'))
		YResponse = make_ResponsePlot(config, lowerData, dictPlot, bosonNameLong, 'Boson_Pt',  10,200,10,0,0,'PF Met','p_T^H','GeV')
		XRange, YVariance = make_ResolutionPlot(config, lowerData, dictPlot, bosonNameLong, 'Boson_Pt',  10,200,10,0,0,'PF Met','p_T^H','GeV')


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
		YResponse = make_ResponsePlot(config, upperData, dictPlot, bosonNameLong, 'Boson_Pt',  10,200,10,0,0,'PF Met','p_T^H','GeV')
		XRange, YVariance = make_ResolutionPlot(config, upperData, dictPlot, bosonNameLong, 'Boson_Pt',  10,200,10,0,0,'PF Met','p_T^H','GeV')

	if 'recoilslimmedMETs_LongZ' in dictPlot:
		bosonNameLong = 'recoilslimmedMETs_LongZ'
		bosonNamePerp = 'recoilslimmedMETs_PerpZ'
		if not os.path.exists((config['outputDir'] + 'ControlPlots/SingleDistributions/')):
			os.makedirs((config['outputDir'] + 'ControlPlots/SingleDistributions/'))
		YResponse = make_ResponsePlot(config, upperData, dictPlot, bosonNameLong, 'Boson_Pt',  10,200,10,0,0,'PF Met','p_T^H','GeV')
		XRange, YVariance = make_ResolutionPlot(config, upperData, dictPlot, bosonNameLong, 'Boson_Pt',  10,200,10,0,0,'PF Met','p_T^H','GeV')

def make_MoreBDTPlots(config, plotData, dictPlot):

	num_bins = 50

	if not os.path.exists((config['outputDir'] + 'CustomPlots/')):
		os.makedirs((config['outputDir'] + 'CustomPlots/'))


	plt.clf()
	plt.hist(plotData[:,dictPlot['Boson_Pt']], num_bins, histtype = 'step', range=[0,200], facecolor='green')
	plt.xlabel(r'$p_T^H\ \mathrm{in\ GeV}$',fontsize = 22)
	plt.ylabel(r'$\mathrm{Events}$',fontsize = 22)
	plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
	plt.plot((5,5),(plt.ylim()[0], plt.ylim()[1]), 'k--')
	plt.plot((15,15),(plt.ylim()[0], plt.ylim()[1]), 'k--')
	plt.plot((25,25),(plt.ylim()[0], plt.ylim()[1]), 'k--')
	plt.plot((35,35),(plt.ylim()[0], plt.ylim()[1]), 'k--')
	plt.plot((45,45),(plt.ylim()[0], plt.ylim()[1]), 'k--')
	plt.plot((60,60),(plt.ylim()[0], plt.ylim()[1]), 'k--')
	plt.plot((80,80),(plt.ylim()[0], plt.ylim()[1]), 'k--')
	plt.plot((100,100),(plt.ylim()[0], plt.ylim()[1]), 'k--')
	plt.plot((150,150),(plt.ylim()[0], plt.ylim()[1]), 'k--')
	plt.plot((200,200),(plt.ylim()[0], plt.ylim()[1]), 'k--')

	plt.tight_layout()
	plt.savefig(config['outputDir'] + "/CustomPlots/PtDistribution.pdf")

	num_bins = 40
	plt.clf()
	plt.hist(plotData[plotData[:,dictPlot['NVertex']]<40,dictPlot['NVertex']], num_bins, histtype = 'step', range=[0,40], facecolor='green')
	plt.xlabel(r'#$PV$',fontsize = 20)
	plt.ylabel(r'$\mathrm{Events}$',fontsize = 22)
	plt.plot((0,0),(plt.ylim()[0], plt.ylim()[1]), 'k--')
	plt.plot((10,10),(plt.ylim()[0], plt.ylim()[1]), 'k--')
	plt.plot((15,15),(plt.ylim()[0], plt.ylim()[1]), 'k--')
	plt.plot((20,20),(plt.ylim()[0], plt.ylim()[1]), 'k--')
	plt.plot((25,25),(plt.ylim()[0], plt.ylim()[1]), 'k--')
	plt.plot((30,30),(plt.ylim()[0], plt.ylim()[1]), 'k--')
	plt.plot((40,40),(plt.ylim()[0], plt.ylim()[1]), 'k--')
	plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
	plt.tight_layout()
	plt.savefig(config['outputDir'] + "/CustomPlots/NPUDistribution.pdf")

	num_bins = 50

	plt.clf()
	index = 0
	for variableName in dictPlot:
		if re.search("sumEtFract", variableName):
			index += 1
			data = plotData[(plotData[:,dictPlot[variableName]]>0) & (plotData[:,dictPlot[variableName]] < 1),dictPlot[variableName]]

			if re.search("PUCorrectedMET", variableName):
				labelName = r'$\mathrm{PU\ Corrected}\,E_T^{miss}$'
			if re.search("slimmedMETsPuppi", variableName):
				labelName = r'$\mathrm{PUPPI}\,E_T^{miss}$'
			if re.search("pfNoPUMET", variableName):
				labelName = r'$\mathrm{No\ PU}\,E_T^{miss}$'
			if re.search("TrackMET", variableName):
				labelName = r'$\mathrm{Track}\,E_T^{miss}$'
			if re.search("pfPUMET", variableName):
				labelName = r'$\mathrm{PU}\,E_T^{miss}$'
			if re.search("slimmedMETs_", variableName):
				labelName = r'$\mathrm{PF}\,E_T^{miss}$'
			plt.hist(data, num_bins, histtype='step', label=labelName, range=[0,1])
			plt.text(0.18*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],(1.05-0.05*index)*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$\mu ($%s$) = %.3f$'%(labelName[:-13]+'$',plotData[:,dictPlot[variableName]].mean()),color = 'k',fontsize=16)
	plt.legend(loc='best', numpoints=1, fontsize = 12)
	plt.xlabel(r'$\sum E_T / \sum E_T^{PF}$',fontsize = 16)
	plt.ylabel(r'$\mathrm{Events}$',fontsize = 18)
	plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
	plt.tight_layout()
	plt.savefig(config['outputDir'] + "/CustomPlots/SumETFraction.png")


	if 'recoilslimmedMETs_Phi' in dictPlot and 'Boson_Phi' in dictPlot and 'PhiCorrectedRecoil_Phi' in dictPlot:
		DPhiPFBoson = plotData[:,dictPlot["Boson_Phi"]]-plotData[:,dictPlot["recoilslimmedMETs_Phi"]]+np.pi - 2.*np.pi*((plotData[:,dictPlot["Boson_Phi"]]-plotData[:,dictPlot["recoilslimmedMETs_Phi"]])>0)
		DPhiMVABoson = plotData[:,dictPlot["Boson_Phi"]]-plotData[:,dictPlot["PhiCorrectedRecoil_Phi"]]+np.pi - 2.*np.pi*((plotData[:,dictPlot["Boson_Phi"]]-plotData[:,dictPlot["PhiCorrectedRecoil_Phi"]])>0)

		plt.clf()
		n, bins, patches = plt.hist(DPhiPFBoson, num_bins, facecolor='green', label=(r'$\Delta(\phi_{PF},\phi_{Z})$'))
		plt.legend(loc='best', numpoints=1)
		plt.xlabel(r'$\phi$',fontsize = 20)
		plt.ylabel('Frequency distribution',fontsize = 20)
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "/CustomPlots/DPhiPFBoson.pdf")

		plt.clf()
		n, bins, patches = plt.hist(DPhiMVABoson, num_bins, facecolor='green', label=(r'$\Delta(\phi_{MVA},\phi_{Z})$'))

		plt.legend(loc='best', numpoints=1)
		plt.xlabel(r'$\phi$',fontsize = 20)
		plt.ylabel('Frequency distribution',fontsize = 20)
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "/CustomPlots/DPhiMVABoson.pdf")

		plt.clf()

		n, bins, patches = plt.hist([DPhiMVABoson,DPhiPFBoson], num_bins, label=[r'$\Delta(\phi_{MVA},\phi_{Z})$',r'$\Delta(\phi_{PF},\phi_{Z})$'])

		plt.legend(loc='best', numpoints=1)
		plt.xlabel(r'$\phi$',fontsize = 20)
		plt.ylabel('Frequency distribution',fontsize = 20)
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "/CustomPlots/DPhiBothBoson.pdf")

		plt.clf()


		slicedDataLow = plotData[60>plotData[:,dictPlot['Boson_Pt']],:]
		slicedDataHigh = plotData[60<=plotData[:,dictPlot['Boson_Pt']],:]
		DPhiPFBosonLow = slicedDataLow[:,dictPlot["Boson_Phi"]]-slicedDataLow[:,dictPlot["recoilslimmedMETs_Phi"]]+np.pi - 2.*np.pi*((slicedDataLow[:,dictPlot["Boson_Phi"]]-slicedDataLow[:,dictPlot["recoilslimmedMETs_Phi"]])>0)
		DPhiPFBosonHigh = slicedDataHigh[:,dictPlot["Boson_Phi"]]-slicedDataHigh[:,dictPlot["recoilslimmedMETs_Phi"]]+np.pi - 2.*np.pi*((slicedDataHigh[:,dictPlot["Boson_Phi"]]-slicedDataHigh[:,dictPlot["recoilslimmedMETs_Phi"]])>0)
		n, bins, patches = plt.hist(DPhiPFBosonLow, num_bins, histtype = 'step', label=(r'$\mathrm{Target}\ T_2 = \phi_Z - \phi_U - \pi$'), range=[-np.pi,np.pi])
		plt.xlabel(r'$\mathrm{Target}\ T_1$',fontsize = 20)
		plt.plot((0.5,0.5),(plt.ylim()[0], plt.ylim()[1]), 'r--', label=r'$\mathrm{Target\ cut} \pm 0.5$')
		plt.plot((-0.5,-0.5),(plt.ylim()[0], plt.ylim()[1]), 'r--')
		plt.ylabel(r'$\mathrm{Events}$',fontsize = 20)
		plt.legend(loc='best', numpoints=1,fontsize = 20)
		plt.text(0.65*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.6*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$p_T^H < 60\,\mathrm{GeV}$',color = 'k',fontsize=16, weight='bold', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
		plt.xlim([-np.pi,np.pi])
		plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "/CustomPlots/TargetPhiLow.pdf")

		plt.clf()

		n, bins, patches = plt.hist(DPhiPFBosonHigh, num_bins, histtype = 'step', label=(r'$\mathrm{Target}\ T_2 =\phi_Z - \phi_U - \pi$'), range=[-np.pi,np.pi])
		plt.xlabel(r'$\mathrm{Target}\ T_1$',fontsize = 20)
		plt.plot((0.5,0.5),(plt.ylim()[0], plt.ylim()[1]), 'r--', label=r'$\mathrm{Target\ cut} \pm 0.5$')
		plt.plot((-0.5,-0.5),(plt.ylim()[0], plt.ylim()[1]), 'r--')
		plt.ylabel(r'$\mathrm{Events}$',fontsize = 20)
		plt.legend(loc='best', numpoints=1, fontsize = 20)
		plt.text(0.65*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.6*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$p_T^H > 60\,\mathrm{GeV}$',color = 'k',fontsize=16, weight='bold', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
		plt.xlim([-np.pi,np.pi])
		plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "/CustomPlots/TargetPhiHigh.pdf")

		plt.clf()

	if 'PhiCorrectedRecoil_LongZ' in dictPlot and 'Boson_Pt' in dictPlot:
		lowerPtCut = 80
		upperPtCut = 100

		plotDataCut = plotData[(lowerPtCut<plotData[:,dictPlot['Boson_Pt']]) & (upperPtCut>plotData[:,dictPlot['Boson_Pt']])]
		BosonPtDistri = plotDataCut[:,dictPlot["Boson_Pt"]]
		PhiPtDistri = -plotDataCut[:,dictPlot["PhiCorrectedRecoil_LongZ"]]
		TargetDistri = BosonPtDistri/PhiPtDistri

		plt.clf()
		n, bins, patches = plt.hist(BosonPtDistri, num_bins, facecolor='green', label=(r'$p_T^H$'))
		plt.legend(loc='best', numpoints=1)
		plt.xlabel(r'$p_t$ in GeV',fontsize = 20)
		plt.ylabel('Frequency distribution',fontsize = 20)
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "/CustomPlots/BosonPt_%i_to_%i.pdf"%(lowerPtCut,upperPtCut))

		plt.clf()
		n, bins, patches = plt.hist([-plotDataCut[:,dictPlot["LongZCorrectedRecoil_LongZ"]],-plotDataCut[:,dictPlot["recoilslimmedMETs_LongZ"]]], num_bins, facecolor='green', label=(r'$\mathrm{MVA Met}$',r'$\mathrm{PF Met}$'), histtype = 'step', range=[lowerPtCut-55,upperPtCut+55])
		plt.legend(loc='best', numpoints=1)
		plt.xlabel(r'$p_t$ in GeV',fontsize = 20)
		plt.ylabel('Frequency distribution',fontsize = 20)
		plt.title(r'$p_T^U$ distribution of METs', fontsize=20)
		plt.text(0.70*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.50*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$\mu_{PF} = %.3f$''\n'r'$\sigma_{PF} = %.3f$''\n'r'$\mu_{MVA} = %.3f$''\n'r'$\sigma_{MVA} = %.3f$''\n'%(-plotDataCut[:,dictPlot["recoilslimmedMETs_LongZ"]].mean(),plotDataCut[:,dictPlot["recoilslimmedMETs_LongZ"]].std(),-plotDataCut[:,dictPlot["LongZCorrectedRecoil_LongZ"]].mean(),plotDataCut[:,dictPlot["LongZCorrectedRecoil_LongZ"]].std()),color = 'k',fontsize=16)
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "/CustomPlots/MetsPtPt_%i_to_%i.pdf"%(lowerPtCut,upperPtCut))


		plt.clf()
		n, bins, patches = plt.hist(PhiPtDistri, num_bins, facecolor='green', label=(r'$\mathrm{MVA Met}_{\phi}$'), range=[lowerPtCut-45,upperPtCut+45])
		plt.legend(loc='best', numpoints=1)
		plt.xlabel(r'$p_t$ in GeV',fontsize = 20)
		plt.ylabel('Frequency distribution',fontsize = 20)
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "/CustomPlots/PhiCorrectedPt_%i_to_%i.pdf"%(lowerPtCut,upperPtCut))

		plt.clf()
		num_bins = 60
		n, bins, patches = plt.hist([BosonPtDistri,PhiPtDistri], num_bins, histtype = 'step', label=[r'$p_T^H$',r'$\mathrm{MVA Met}_{\phi}$'], range=[lowerPtCut-50,upperPtCut+50])

		plt.legend(loc='best', numpoints=1)
		plt.xlabel(r'$p_T$ in GeV',fontsize = 20)
		plt.ylabel('Frequency distribution',fontsize = 20)
		plt.text(0.65*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.4*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$(80 < p_T^H < 100)\,\mathrm{GeV}$',color = 'k',fontsize=16, weight='bold', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "/CustomPlots/BosonAndPhiPt_%i_to%i.pdf"%(lowerPtCut,upperPtCut))

		plt.clf()

		num_bins = 50
		n, bins, patches = plt.hist(TargetDistri, num_bins, histtype = 'step', label=(r'$\mathrm{Target}\ T_2 = -p_T^H / U_{||,\phi}^{\mathrm{MVA}}$'), range=[0.5,2])
		plt.xlabel(r'$\mathrm{Target}\ T_2$',fontsize = 18)
		plt.plot((0.7,0.7),(plt.ylim()[0], plt.ylim()[1]), 'k--', label=r'$\mathrm{Target\ Cut}\ 1 \pm 0.3$')
		plt.plot((1.3,1.3),(plt.ylim()[0], plt.ylim()[1]), 'k--')
		plt.plot((0.6,0.6),(plt.ylim()[0], plt.ylim()[1]), 'g--', label=r'$\mathrm{Target\ Cut}\ 1 \pm 0.4$')
		plt.plot((1.4,1.4),(plt.ylim()[0], plt.ylim()[1]), 'g--')
		plt.plot((0.5,0.5),(plt.ylim()[0], plt.ylim()[1]), 'b--', label=r'$\mathrm{Target\ Cut}\ 1 \pm 0.5$')
		plt.plot((1.5,1.5),(plt.ylim()[0], plt.ylim()[1]), 'b--')
		plt.plot((0.3,0.3),(plt.ylim()[0], plt.ylim()[1]), 'r--', label=r'$\mathrm{Target\ Cut}\ 1 \pm 0.7$')
		plt.plot((1.7,1.7),(plt.ylim()[0], plt.ylim()[1]), 'r--')
		plt.ylabel(r'$\mathrm{Events}$',fontsize = 18)
		plt.legend(loc='upper right', numpoints=1, fontsize = 14)
		plt.text(0.65*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.4*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$(80 < p_T^H < 100)\,\mathrm{GeV}$',color = 'k',fontsize=16, weight='bold', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
		plt.text(0.35*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.15*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$\mu_{empiric} = %.3f$''\n'%TargetDistri.mean(),color = 'k',fontsize=16)
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "/CustomPlots/Target_%i_to_%i.pdf"%(lowerPtCut,upperPtCut))


		plt.clf()

	if 'recoilslimmedMETs_LongZ' in dictPlot and 'Boson_Pt' in dictPlot and 'LongZCorrectedRecoil_LongZ' in dictPlot and 'LongZCorrectedRecoil_PerpZ' in dictPlot:
		plt.clf()
		MVARecoilDiffDistri = np.sqrt((plotData[:,dictPlot["LongZCorrectedRecoil_LongZ"]]+plotData[:,dictPlot["Boson_Pt"]])**2+(plotData[:,dictPlot["LongZCorrectedRecoil_PerpZ"]])**2)
		PFRecoilDiffDistri = np.sqrt((plotData[:,dictPlot["recoilslimmedMETs_LongZ"]]+plotData[:,dictPlot["Boson_Pt"]])**2 +(plotData[:,dictPlot["recoilslimmedMETs_PerpZ"]])**2)
		n, bins, patches = plt.hist([PFRecoilDiffDistri,MVARecoilDiffDistri], num_bins, label=[r'$Reco_{Recoil}$',r'$MVA_{Recoil}$'], range=[0,80])
		plt.xlabel(r"MET in GeV",fontsize = 20)
		plt.ylabel("Frequency Distribution", fontsize = 20)
		plt.legend(loc='best', numpoints=1, fontsize = 20)
		plt.text(0.70*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.50*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$\mu_{Reco} = %.3f$''\n'r'$\mu_{MVA} = %.3f$''\n'%(PFRecoilDiffDistri.mean(),MVARecoilDiffDistri.mean()),color = 'k',fontsize=20)
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "/CustomPlots/Diff_Recoil_BosonPt.pdf")
		plt.clf()

	if 'recoilslimmedMETs_Phi' in dictPlot and 'Boson_Phi' in dictPlot and 'LongZCorrectedRecoil_Phi' in dictPlot:
		slicedDataLowest = plotData[(200>plotData[:,dictPlot['Boson_Pt']]) &(60<plotData[:,dictPlot['Boson_Pt']]),:]
		DPhiPFBosonLowest = slicedDataLowest[:,dictPlot["Boson_Phi"]]-slicedDataLowest[:,dictPlot["recoilslimmedMETs_Phi"]]+np.pi - 2.*np.pi*((slicedDataLowest[:,dictPlot["Boson_Phi"]]-slicedDataLowest[:,dictPlot["recoilslimmedMETs_Phi"]])>0) + 2.*np.pi*((slicedDataLowest[:,dictPlot["Boson_Phi"]]-slicedDataLowest[:,dictPlot["recoilslimmedMETs_Phi"]])<-2.*np.pi)

		DPhiMVABosonLowest = slicedDataLowest[:,dictPlot["Boson_Phi"]]-slicedDataLowest[:,dictPlot["LongZCorrectedRecoil_Phi"]]+np.pi - 2.*np.pi*((slicedDataLowest[:,dictPlot["Boson_Phi"]]-slicedDataLowest[:,dictPlot["LongZCorrectedRecoil_Phi"]])>0) + 2.*np.pi*((slicedDataLowest[:,dictPlot["Boson_Phi"]]-slicedDataLowest[:,dictPlot["LongZCorrectedRecoil_Phi"]])<-2.*np.pi) 
		#n, bins, patches = plt.hist([DPhiPFBosonLowest,DPhiMVABosonLowest], num_bins, histtype = 'step', label=[r'\Delta(PF) = $\phi_Z - \phi_U^{PF} - \pi$',r'\Delta(MVA) = $\phi_Z - \phi_U^{MVA} - \pi$'], range=[-np.pi,np.pi])
		n, bins, patches = plt.hist([DPhiPFBosonLowest,DPhiMVABosonLowest], num_bins, histtype = 'step', label=[r'\Delta(PF) = $\phi_Z - \phi_U^{PF} - \pi$',r'\Delta(MVA) = $\phi_Z - \phi_U^{MVA} - \pi$'])
		plt.xlabel(r'$\phi_Z - \phi_U - \pi$',fontsize = 20)
		plt.ylabel('Frequency distribution',fontsize = 20)
		plt.legend(loc='best', numpoints=1)
		plt.text(0.65*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.6*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$(15<p_T^H < 25)\,\mathrm{GeV}$',color = 'k',fontsize=16, weight='bold', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
		#plt.xlim([-np.pi,np.pi])
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "/CustomPlots/AngularDifferenceMVAAndPF.pdf")

		plt.clf()

	if 'probChiSquare' in dictPlot and 'probChiSquarePf' in dictPlot:
		probChiMVA = plotData[:,dictPlot["probChiSquare"]]
		probChiPf = plotData[:,dictPlot["probChiSquarePf"]]

		n, bins, patches = plt.hist([probChiMVA,probChiPf], num_bins, histtype = 'step', label=[r'MVA Met',r'PF Met'])
		plt.xlabel(r"$prob(\chi^2)$", fontsize = 20)
		plt.ylabel("Events", fontsize = 20)
		plt.legend(loc='lower right', numpoints=1, fontsize = 20)
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "/CustomPlots/ChiSquareProbDistri.pdf")
		plt.clf()

	if 'Boson_Pt' in dictPlot and 'LongZCorrectedRecoil_LongZ' in dictPlot:
		#centralData = data[((weightedMean(data[:,0],data[:,1])-4*weightedStd(data[:,0],data[:,1]))<data[:,0]) & (data[:,0] <(weightedMean(data[:,0],data[:,1])+4*weightedStd(data[:,0],data[:,1])))]
		AlternativeDistri = plotData[(100<plotData[:,dictPlot["Boson_Pt"]]) & (110>plotData[:,dictPlot["Boson_Pt"]])]
		BosonPtSelected = AlternativeDistri[:,dictPlot["Boson_Pt"]]
		UParalSelected = -AlternativeDistri[:,dictPlot["LongZCorrectedRecoil_LongZ"]]
		plt.clf()
		n, bins, patches = plt.hist([BosonPtSelected,UParalSelected], num_bins, range=[40,200], histtype = 'step', label=[r'$p_T^H$',r'$U_{||}$'])
		plt.xlabel(r"$p_T^H$", fontsize = 20)
		plt.ylabel("Events", fontsize = 20)
		plt.legend(loc='best', numpoints=1, fontsize = 20)
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "/CustomPlots/AsymmetryInput.pdf")
		plt.clf()

		BosonOverU = BosonPtSelected/UParalSelected
		UOverBoson = UParalSelected/BosonPtSelected

		plt.clf
		n, bins, patches = plt.hist([UOverBoson,BosonOverU], num_bins, normed=True,range=[0,2.5], histtype = 'step', label=[r'$U_{||}/p_T^H$',r'$p_T^H/U_{||}$'])
		yUOverB = mlab.normpdf(bins, UOverBoson.mean(), UOverBoson.std())
		yBOverU = mlab.normpdf(bins, BosonOverU.mean(), BosonOverU.std())
		#plt.plot(bins, yUOverB, 'b--')
		#plt.plot(bins, yBOverU, 'g--')
		plt.xlabel(r"$Response$", fontsize = 20)
		plt.ylabel("Distribution function", fontsize = 20)
		plt.legend(loc='best', numpoints=1, fontsize = 20)
		plt.text(0.60*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.30*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$Events=%i$''\n'r'$<p_T^H/U_{||}> = %.3f$''\n'r'$<U_{||}/p_T^H> = %.3f$''\n'%(AlternativeDistri.shape[0],UOverBoson.mean(),BosonOverU.mean()),color = 'k',fontsize=20)
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "/CustomPlots/AsymmetryOutput.pdf")
		plt.clf()


	bashCommand = 'cp index.php %s/CustomPlots/'%(config['outputDir'])
	os.system(bashCommand)

	print("Custom plots created")



def make_PhiVariancePlot(config, plotData, dictPlot, targetvariable, ptmin, ptmax, xlabelname = ''):

	if xlabelname == '':
		xlabelname = targetvariable
	num_bins = 50
	if 'Boson_Phi' in dictPlot:
		histDataPhi = plotData[:,dictPlot['Boson_Phi']] + math.pi - plotData[:,dictPlot[targetvariable]]
		print(targetvariable,' phi shape: ',histDataPhi.shape)
		for event in range(0,histDataPhi.shape[0]):
			if histDataPhi[event] > math.pi:
				histDataPhi[event] -= 2*math.pi
			if histDataPhi[event] < -math.pi:
				histDataPhi[event] += 2*math.pi
		MSE = (histDataPhi**2).mean()
		n, bins, patches = plt.hist(histDataPhi, num_bins, facecolor='green', alpha=0.5)
		plt.xlabel('Variance %s from true Boson Phi (%i < Boson Pt < %i)GeV. MSE: %f'%(xlabelname, ptmin, ptmax, MSE),fontsize = 20)
		plt.ylabel('Entries',fontsize = 20)
		plt.savefig(config['outputDir'] + "/CustomPlots/PhiVariance%s_%ito%iPt.pdf"%(xlabelname,ptmin,ptmax))
		plt.clf()
		print('MSE %s %ito%iGeV: '%(xlabelname,ptmin,ptmax),MSE)

		# normal distribution center at x=0 and y=5
		plt.hist2d(plotData[:,dictPlot['Boson_Phi']], histDataPhi,bins = 80, norm=LogNorm())
		#plt.ylim([-0.25,0.25])
		plt.xlabel('Boson Phi (%i < Boson Pt < %i)GeV'%(ptmin, ptmax),fontsize = 20)
		plt.ylabel('Variance of (Prediction-Target) %s'%xlabelname,fontsize = 20)
		plt.colorbar()
		plt.tight_layout()
		plt.savefig(config['outputDir'] + "/CustomPlots/Variance2D_%s%ito%iGeV.pdf"%(xlabelname,ptmin,ptmax))
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

	#Jet study plots, can be enabled in the config
	if plotconfig['plotJetStudyPlots'] and 'Jet0_Pt' in dictPlot and 'Jet1_Pt' in dictPlot:
		make_JetStudyPlots(config, plotData, dictPlot)

	if 'skipControlPlots' in config:
		if config['skipControlPlots']:
			return True
	#Boson Pt
	#comparisonMinus.xlabel('Boson_Pt')
	#comparisonOver.xlabel('Boson_Pt')
	#comparisonMinus.ylabel('|Boson_Pt-Prediction|')
	#comparisonOver.ylabel('Prediction/Boson_Pt')

	bosonmin = [0,0,plotconfig['BosonCut']]
	bosonmax = [plotData[:,dictPlot['Boson_Pt']].max(),plotconfig['BosonCut'],plotData[:,dictPlot['Boson_Pt']].max()]

	#BDT Performance Plots
#Phi ControlPlots

	binRangesNPV = np.array([0,10,15,20,25,30,40])
	for i, min in enumerate(bosonmin):
		slicedData = plotData[min<=plotData[:,dictPlot['Boson_Pt']],:]
		slicedData = slicedData[bosonmax[i]>=slicedData[:,dictPlot['Boson_Pt']],:]
		if plotconfig['plotBDTPerformance']:
			if plotconfig['plotAdditionalBDTPlots']:
				make_PhiVariancePlot(config, slicedData,dictPlot,'PhiCorrectedRecoil_Phi', min, bosonmax[i], 'BDT Phi')


				#BDT

			if 'LongZCorrectedRecoil_LongZ' in dictPlot:
				if 'inputLabel' in config:
					make_ControlPlots(config, slicedData, dictPlot, 'LongZCorrectedRecoil', 'NVertex', 5,40,5,min,bosonmax[i],config['inputLabel'][0],'\# PV','',binRanges=binRangesNPV)
				else:
					make_ControlPlots(config, slicedData, dictPlot, 'LongZCorrectedRecoil', 'NVertex', 5,40,5,min,bosonmax[i],r'$\mathrm{MVA}\ E_T^{miss}$','\# PV','',binRanges=binRangesNPV)


			for index in range(len(config['inputFile'])-1):	
				if 'inputLabel' in config:
					make_ControlPlots(config, slicedData, dictPlot, 'V%iLongZCorrectedRecoil'%index, 'NVertex', 5,40,5,min,bosonmax[i],config['inputLabel'][index+1],'\# PV','',binRanges=binRangesNPV)
				else:
					make_ControlPlots(config, slicedData, dictPlot, 'V%iLongZCorrectedRecoil'%index, 'NVertex', 5,40,5,min,bosonmax[i],'MVA Met %i'%(index+2),'\# PV','',binRanges=binRangesNPV)


			if plotconfig['plotPhiCorrected']:
				if 'PhiCorrectedRecoil_LongZ' in dictPlot:
					if 'inputLabel' in config:
						make_ControlPlots(config, slicedData, dictPlot, 'PhiCorrectedRecoil', 'NVertex', 5,40,5,min,bosonmax[i],config['inputLabel'][0],'\# PV','',binRanges=binRangesNPV)
					else:
						make_ControlPlots(config, slicedData, dictPlot, 'PhiCorrectedRecoil', 'NVertex', 5,40,5,min,bosonmax[i],r'$\phi - \mathrm{MVA}\ E_T^{miss}$','\# PV','',binRanges=binRangesNPV)


				for index in range(len(config['inputFile'])-1):
					if 'inputLabel' in config:
						make_ControlPlots(config, slicedData, dictPlot, 'V%iPhiCorrectedRecoil'%index, 'NVertex', 5,40,5,min,bosonmax[i],config['inputLabel'][index+1],'\# PV','',binRanges=binRangesNPV)
					else:
						make_ControlPlots(config, slicedData, dictPlot, 'V%iPhiCorrectedRecoil'%index, 'NVertex', 5,40,5,min,bosonmax[i],'MVA Met %i'%(index+2),'\# PV','',binRanges=binRangesNPV)


		#slimmedMet
		if 'recoilslimmedMETs_LongZ' in dictPlot:
			make_ControlPlots(config, slicedData, dictPlot, 'recoilslimmedMETs', 'NVertex', 5,40,5,min,bosonmax[i],r'$\mathrm{PF}\ E_T^{miss}$','\# PV','',binRanges=binRangesNPV)


		#Puppi-Met
		if plotconfig['plotPuppiPerformance']:
			if 'recoilslimmedMETsPuppi_LongZ' in dictPlot:
				make_ControlPlots(config, slicedData, dictPlot, 'recoilslimmedMETsPuppi', 'NVertex', 5,40,5,min,bosonmax[i],r'$\mathrm{PUPPI}\ E_T^{miss}$','\# PV','',binRanges=binRangesNPV)


		if plotconfig['plotAdditionalMETPerformance']:
			if 'recoilpatpfTrackMET_LongZ' in dictPlot:
				make_ControlPlots(config, slicedData, dictPlot, 'recoilpatpfTrackMET', 'NVertex', 5,40,5,min,bosonmax[i],r'$\mathrm{Track}\ E_T^{miss}$','\# PV','',binRanges=binRangesNPV)
			if 'recoilpatpfPUMET_LongZ' in dictPlot:
				make_ControlPlots(config, slicedData, dictPlot, 'recoilpatpfPUMET', 'NVertex', 5,40,5,min,bosonmax[i],r'$\mathrm{PU}\ E_T^{miss}$','\# PV','',binRanges=binRangesNPV)
			if 'recoilpatpfPUCorrectedMET_LongZ' in dictPlot:
				make_ControlPlots(config, slicedData, dictPlot, 'recoilpatpfPUCorrectedMET', 'NVertex', 5,40,5,min,bosonmax[i],r'$\mathrm{PU Corrected}\ E_T^{miss}$','\# PV','',binRanges=binRangesNPV)
			if 'recoilpatpfNoPUMET_LongZ' in dictPlot:
				make_ControlPlots(config, slicedData, dictPlot, 'recoilpatpfNoPUMET', 'NVertex', 5,40,5,min,bosonmax[i],r'$\mathrm{No PU}\ E_T^{miss}$','\# PV','',binRanges=binRangesNPV)

		plt.clf()
		plt.figure(4)
		#plt.xlabel(r'#PV ($%i\,\mathrm{GeV} < p_T^H < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]),fontsize = 20)
		if config['fixedRange']:
			plt.ylim([0.6,1.05])
		if min > 0:
			plt.text(0.60*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.15*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$p_T^H > %i\,\mathrm{GeV}$'%min,color = 'k',fontsize=24, weight='bold', bbox=dict(facecolor='none', edgecolor='black', pad=10.0))
		plt.xlabel(r'#PV',fontsize = 20)
		plt.ylabel(r'$\left<U_{||} / -p_T^H\right>$',fontsize = 22)
		#plt.ylabel(r'$<-p_T^H/U_{\|}>$',fontsize = 20)
		plt.title(r'$\mathrm{Response}\ U_{||}$', fontsize = 28, y = 1.02)
		legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
		plt.plot((plt.xlim()[0], plt.xlim()[1]), (1, 1), 'k--')
		plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
		plt.tight_layout()
		plt.savefig(config['outputDir'] + 'Response_%ito%iGeV_vs_NVertex.pdf'%(min,bosonmax[i]))
		plt.clf()

		plt.figure(5)
		#plt.xlabel(r'#PU ($%i\,\mathrm{GeV} < p_T^H < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]),fontsize = 20)
		if config['fixedRange']:
			plt.ylim([12,30])
		if min > 0:
			plt.text(0.20*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.10*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$p_T^H > %i\,\mathrm{GeV}$'%min,color = 'k',fontsize=24, weight='bold', bbox=dict(facecolor='none', edgecolor='black', pad=10.0))
		plt.xlabel(r'#PV',fontsize = 20)
		plt.ylabel(r'$\sigma(U_{||} + p_T^H)\ \mathrm{in\ GeV}$',fontsize = 22)
		plt.title(r'$\mathrm{Resolution}\ U_{||}$', fontsize = 28, y = 1.02)
		legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
		plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
		plt.tight_layout()
		plt.savefig(config['outputDir'] + 'Resolution_%ito%iGeV_vs_NVertex.pdf'%(min,bosonmax[i]))
		plt.clf()

		plt.figure(6)
		if config['fixedRange']:
			plt.ylim([12,30])
		if min > 0:
			plt.text(0.20*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.10*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$p_T^H > %i\,\mathrm{GeV}$'%min,color = 'k',fontsize=24, weight='bold', bbox=dict(facecolor='none', edgecolor='black', pad=10.0))
		plt.xlabel(r'#PV',fontsize = 20)
		plt.ylabel(r'$\mathrm{Resolution} / \mathrm{Response}$',fontsize = 22)
		plt.title(r'$\mathrm{Response\ Corrected}$', fontsize = 28, y = 1.02)
		legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
		plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
		plt.tight_layout()
		plt.savefig(config['outputDir'] + 'ResponseCorrected_%ito%iGeV_vs_NVertex.pdf'%(min,bosonmax[i]))
		plt.clf()

		plt.figure(7)
		if config['fixedRange']:
			plt.ylim([-1,1])
		if min > 0:
			plt.text(0.60*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.15*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$p_T^H > %i\,\mathrm{GeV}$'%min,color = 'k',fontsize=24, weight='bold', bbox=dict(facecolor='none', edgecolor='black', pad=10.0))
		legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
		plt.xlabel(r'#PV',fontsize = 20)
		plt.ylabel(r'$\left<U_\bot\right>$',fontsize = 22)
		plt.title(r'$\mathrm{Response}\ U_\bot$', fontsize = 28, y = 1.02)
		plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
		plt.tight_layout()
		plt.savefig(config['outputDir'] + 'ResponsePerp_%ito%iGeV_vs_NVertex.pdf'%(min,bosonmax[i]))
		plt.clf()

		plt.figure(8)
		#plt.xlabel(r'#PU ($%i\,\mathrm{GeV} < p_T^H < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]),fontsize = 20)
		if config['fixedRange']:
			plt.ylim([12,30])
		if min > 0:
			plt.text(0.60*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.15*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$p_T^H > %i\,\mathrm{GeV}$'%min,color = 'k',fontsize=24, weight='bold', bbox=dict(facecolor='none', edgecolor='black', pad=10.0))
		plt.xlabel(r'#PV',fontsize = 20)
		plt.ylabel(r'$\sigma(U_\bot)\ \mathrm{in\ GeV}$',fontsize = 22)
		plt.title(r'$\mathrm{Resolution}\ U_\bot$', fontsize = 28, y = 1.02)
		legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
		plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
		plt.tight_layout()
		plt.savefig(config['outputDir'] + 'ResolutionPerp_%ito%iGeV_vs_NVertex.pdf'%(min,bosonmax[i]))
		plt.clf()

		plt.figure(9)
		if config['fixedRange']:
			plt.ylim([0.8,2.3])
		if min > 0:
			plt.text(0.60*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.15*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$p_T^H > %i\,\mathrm{GeV}$'%min,color = 'k',fontsize=24, weight='bold', bbox=dict(facecolor='none', edgecolor='black', pad=10.0))
		#plt.xlabel(r'#PU ($%i\,\mathrm{GeV} < p_T^H < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]),fontsize = 20)
		plt.xlabel(r'#PV',fontsize = 20)
		plt.ylabel(r'$\left<E_{t,||}^{miss}/E_{t,gen}^{miss}\right>$',fontsize = 22)
		plt.title(r'$\mathrm{Response}\ E_{T,||}^{miss}$', fontsize = 28, y = 1.02)
		plt.plot((plt.xlim()[0], plt.xlim()[1]), (1, 1), 'k--')
		legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
		plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
		plt.tight_layout()
		plt.savefig(config['outputDir'] + 'METResponse_%ito%iGeV_vs_NVertex.pdf'%(min,bosonmax[i]))
		plt.clf()

		plt.figure(10)
		if config['fixedRange']:
			plt.ylim([12.5,32.5])
		if min > 0:
			plt.text(0.60*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.15*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$p_T^H > %i\,\mathrm{GeV}$'%min,color = 'k',fontsize=24, weight='bold', bbox=dict(facecolor='none', edgecolor='black', pad=10.0))
		plt.ylabel(r'$\sigma(E_{t,||}^{miss}-E_{t,gen}^{miss})\ \mathrm{in\ GeV}$',fontsize = 22)
		#plt.xlabel(r'#PU ($%i\,\mathrm{GeV} < p_T^H < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]),fontsize = 20)
		#location lower right
		plt.xlabel(r'#PV',fontsize = 20)

		plt.title(r'$\mathrm{Resolution}\ E_{T,||}^{miss}$', fontsize = 28, y = 1.02)
		legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
		plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
		plt.tight_layout()
		plt.savefig(config['outputDir'] + 'METResolution_%ito%iGeV_vs_NVertex.pdf'%(min,bosonmax[i]))
		plt.clf()
		plt.figure(0)

		plt.figure(11)
		if config['fixedRange']:
			plt.ylim([12.5,32.5])
		if min > 0:
			plt.text(0.60*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.15*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$p_T^H > %i\,\mathrm{GeV}$'%min,color = 'k',fontsize=24, weight='bold', bbox=dict(facecolor='none', edgecolor='black', pad=10.0))
		plt.ylabel(r'$\sigma(E_{t\bot}^{miss})\ \mathrm{in\ GeV}$',fontsize = 22)
		plt.xlabel(r'#PV',fontsize = 20)
		plt.title(r'$\mathrm{Resolution}\ E_{T,\bot}^{miss}$', fontsize = 28, y = 1.02)
		legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
		plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
		plt.tight_layout()
		plt.savefig(config['outputDir'] + 'METResolutionPerp_%ito%iGeV_vs_NVertex.pdf'%(min,bosonmax[i]))
		plt.clf()

		plt.clf()
		plt.figure(12)
		if config['fixedRange']:
			plt.ylim([0.0,1.3])
		if min > 0:
			plt.text(0.60*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.15*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$p_T^H > %i\,\mathrm{GeV}$'%min,color = 'k',fontsize=24, weight='bold', bbox=dict(facecolor='none', edgecolor='black', pad=10.0))
		plt.xlabel(r'#PV',fontsize = 20)
		#plt.ylabel(r'$<U_{\|} / -p_T^H>$',fontsize = 20)
		plt.ylabel(r'$\left<-p_T^H/U_{||}\right>$',fontsize = 22)
		plt.title(r'$\mathrm{Response^{-1}}\ U_{||}$', fontsize = 28, y = 1.02)
		legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
		plt.plot((plt.xlim()[0], plt.xlim()[1]), (1, 1), 'k--')
		plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
		plt.tight_layout()
		plt.savefig(config['outputDir'] + 'ResponseInverted_%ito%iGeV_vs_NVertex.pdf'%(min,bosonmax[i]))
		plt.clf()

		plt.figure(13)
		if config['fixedRange']:
			plt.ylim([0,2])
		if min > 0:
			plt.text(0.60*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],0.15*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],r'$p_T^H > %i\,\mathrm{GeV}$'%min,color = 'k',fontsize=24, weight='bold', bbox=dict(facecolor='none', edgecolor='black', pad=10.0))
		#plt.xlabel(r'#PU ($%i\,\mathrm{GeV} < p_T^H < %i\,\mathrm{GeV}$)'%(min,bosonmax[i]),fontsize = 20)
		plt.xlabel(r'#PV',fontsize = 20)
		plt.ylabel(r'$\left<E_{t,gen}^{miss}/E_{t,||}^{miss}\right>$',fontsize = 22)
		plt.title(r'$\mathrm{Response^{-1}}\ E_{T,||}^{miss}$', fontsize = 28, y = 1.02)
		plt.plot((plt.xlim()[0], plt.xlim()[1]), (1, 1), 'k--')
		legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
		plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
		plt.tight_layout()
		plt.savefig(config['outputDir'] + 'METResponseInverted_%ito%iGeV_vs_NVertex.pdf'%(min,bosonmax[i]))
		plt.clf()
		plt.figure(0)
#additional BDTPlots
		if plotconfig['plotAdditionalBDTPlots']:
			make_PhiVariancePlot(config, slicedData,dictPlot,'recoilslimmedMETs_Phi', min, bosonmax[i], 'PF Phi')

			make_PhiVariancePlot(config, slicedData,dictPlot,'recoilslimmedMETsPuppi_Phi', min, bosonmax[i],'PUPPI Phi')



	#Boson PT


	binRangesPt = np.array([5,15,25,35,45,60,80,100,150,200])
	#BDT Performance
	if plotconfig['plotBDTPerformance']:
#MVAMet
		if 'LongZCorrectedRecoil_LongZ' in dictPlot:
			#If labeling is activated
			if 'inputLabel' in config:
				make_ControlPlots(config, plotData, dictPlot, 'LongZCorrectedRecoil', 'Boson_Pt', 10,200,10,0,0,config['inputLabel'][0],'p_T^H','GeV', binRanges=binRangesPt)
			else:
				make_ControlPlots(config, plotData, dictPlot, 'LongZCorrectedRecoil', 'Boson_Pt', 10,200,10,0,0,r'$\mathrm{MVA}\ E_T^{miss}$','p_T^H','GeV', binRanges=binRangesPt)

			#for additional inputs
			for index in range(len(config['inputFile'])-1):
				if 'inputLabel' in config:
					make_ControlPlots(config, plotData, dictPlot, 'V%iLongZCorrectedRecoil'%index, 'Boson_Pt', 10,200,10,0,0,config['inputLabel'][index+1],'p_T^H','GeV', binRanges=binRangesPt)
				else:
					make_ControlPlots(config, plotData, dictPlot, 'V%iLongZCorrectedRecoil'%index, 'Boson_Pt', 10,200,10,0,0,'MVA Met %i'%(index+2),'p_T^H','GeV', binRanges=binRangesPt)

		#Phi Corrected MVAMET Plots
		if plotconfig['plotPhiCorrected']:
			if 'PhiCorrectedRecoil_LongZ' in dictPlot:
				if 'inputLabel' in config:
					make_ControlPlots(config, plotData, dictPlot, 'PhiCorrectedRecoil', 'Boson_Pt', 10,200,10,0,0,config['inputLabel'][0],'p_T^H','GeV', binRanges=binRangesPt)
				else:
					make_ControlPlots(config, plotData, dictPlot, 'PhiCorrectedRecoil', 'Boson_Pt', 10,200,10,0,0,r'$\phi - \mathrm{MVA}\ E_T^{miss}$','p_T^H','GeV', binRanges=binRangesPt)

				for index in range(len(config['inputFile'])-1):
					if 'inputLabel' in config:
						make_ControlPlots(config, plotData, dictPlot, 'V%iPhiCorrectedRecoil'%index, 'Boson_Pt', 10,200,10,0,0,config['inputLabel'][index+1],'p_T^H','GeV', binRanges=binRangesPt)
					else:
						make_ControlPlots(config, plotData, dictPlot, 'V%iPhiCorrectedRecoil'%index, 'Boson_Pt', 10,200,10,0,0,'MVA Met %i'%(index+2),'p_T^H','GeV', binRanges=binRangesPt)

	#PF Met Control Plots Pt
	if 'recoilslimmedMETs_LongZ' in dictPlot:
		make_ControlPlots(config, plotData, dictPlot, 'recoilslimmedMETs', 'Boson_Pt', 10,200,10,0,0,r'$\mathrm{PF}\ E_T^{miss}$','p_T^H','GeV', binRanges=binRangesPt)


	#Puppi Met Control Plots Pt
	if plotconfig['plotPuppiPerformance']:
		if 'recoilslimmedMETsPuppi_LongZ' in dictPlot:
			make_ControlPlots(config, plotData, dictPlot, 'recoilslimmedMETsPuppi', 'Boson_Pt', 10,200,10,0,0,r'$\mathrm{PUPPI}\ E_T^{miss}$','p_T^H','GeV', binRanges=binRangesPt)


	if plotconfig['plotAdditionalMETPerformance']:
		if 'recoilpatpfTrackMET_LongZ' in dictPlot:
			make_ControlPlots(config, plotData, dictPlot, 'recoilpatpfTrackMET', 'Boson_Pt', 10,200,10,0,0,r'$\mathrm{Track}\ E_T^{miss}$','p_T^H','GeV', binRanges=binRangesPt)
		if 'recoilpatpfPUMET_LongZ' in dictPlot:
			make_ControlPlots(config, plotData, dictPlot, 'recoilpatpfPUMET', 'Boson_Pt', 10,200,10,0,0,r'$\mathrm{PU}\ E_T^{miss}$','p_T^H','GeV', binRanges=binRangesPt)
		if 'recoilpatpfPUCorrectedMET_LongZ' in dictPlot:
			make_ControlPlots(config, plotData, dictPlot, 'recoilpatpfPUCorrectedMET', 'Boson_Pt', 10,200,10,0,0,r'$\mathrm{PU Corrected}\ E_T^{miss}$','p_T^H','GeV', binRanges=binRangesPt)
		if 'recoilpatpfNoPUMET_LongZ' in dictPlot:
			make_ControlPlots(config, plotData, dictPlot, 'recoilpatpfNoPUMET', 'Boson_Pt', 10,200,10,0,0,r'$\mathrm{No PU}\ E_T^{miss}$','p_T^H','GeV', binRanges=binRangesPt)

	plt.figure(4)
	legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
	#legend.get_frame().set_alpha(0.5)
	plt.ylabel(r'$\left<U_{||} / -p_T^H\right>$', fontsize = 22)
	#plt.ylabel(r'$<-p_T^H/U_{\|}>$', fontsize = 20)
	plt.xlabel(r'$p_T^H\ \mathrm{in\ GeV}$', fontsize = 22)
	plt.title(r'$\mathrm{Response}\ U_{||}$', fontsize= 28, y = 1.02)
	plt.plot((plt.xlim()[0], plt.xlim()[1]), (1, 1), 'k--')
	if config['fixedRange']:
		plt.ylim([0.6,1.05])

	plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
	plt.tight_layout()
	plt.savefig(config['outputDir'] + 'Response_vs_BosonPt.pdf')
	plt.clf()
	plt.figure(5)
	plt.xlabel(r'$p_T^H\ \mathrm{in\ GeV}$', fontsize = 22)
	plt.ylabel(r'$\sigma(U_{||} + p_T^H)\ \mathrm{in\ GeV}$', fontsize = 22)
	plt.title(r'$\mathrm{Resolution}\ U_{||}$', fontsize = 28, y = 1.02)
	legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
	#legend.get_frame().set_alpha(0.5)
	if config['fixedRange']:
		plt.ylim([12,30])
	plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
	plt.tight_layout()
	plt.savefig(config['outputDir'] + 'Resolution_vs_BosonPt.pdf')
	plt.clf()
	plt.figure(6)
	plt.xlabel(r'$p_T^H\ \mathrm{in\ GeV}$', fontsize = 22)
	plt.ylabel(r'$\mathrm{Resolution / Response}$',fontsize = 22)
	plt.title(r'$\mathrm{Response\ Corrected}$', fontsize = 28, y = 1.02)
	legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
	if config['fixedRange']:
		plt.ylim([12,30])
	plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
	#legend.get_frame().set_alpha(0.5)
	plt.tight_layout()
	plt.savefig(config['outputDir'] + 'ResponseCorrected_vs_BosonPt.pdf')
	plt.clf()
	plt.figure(7)
	legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
	#legend.get_frame().set_alpha(0.5)
	plt.xlabel(r'$p_T^H\ \mathrm{in\ GeV}$', fontsize = 22)
	plt.ylabel(r'$\left<U_\bot\right>$',fontsize = 22)
	plt.title(r'$\mathrm{Response}\ U_\bot$', fontsize = 28, y = 1.02)
	if config['fixedRange']:
		plt.ylim([-1,1])
	plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
	plt.tight_layout()
	plt.savefig(config['outputDir'] + 'ResponsePerp_vs_BosonPt.pdf')
	plt.clf()
	plt.figure(8)
	plt.xlabel(r'$p_T^H\ \mathrm{in\ GeV}$', fontsize = 22)
	plt.ylabel(r'$\sigma(U_\bot)\ \mathrm{in\ GeV}$',fontsize = 22)
	plt.title(r'$\mathrm{Resolution}\ U_\bot$', fontsize = 28, y = 1.02)
	legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
	#legend.get_frame().set_alpha(0.5)
	if config['fixedRange']:
		plt.ylim([12,30])
	plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
	plt.tight_layout()
	plt.savefig(config['outputDir'] + 'ResolutionPerp_vs_BosonPt.pdf')
	plt.clf()
	plt.figure(9)
	legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
	#legend.get_frame().set_alpha(0.5)
	plt.xlabel(r'$p_T^H\ \mathrm{in\ GeV}$', fontsize = 22)
	#plt.ylabel('MET_Long/genMET',fontsize = 20)
	plt.ylabel(r'$\left<E_{T,||}^{miss}/E_{T,gen}^{miss}\right>$',fontsize = 22)
	#plt.ylabel(r'$\ensuremath{{\not\mathrel{E}}_T}$',fontsize = 20)
	plt.title(r'$\mathrm{Response}\ E_{T,||}^{miss}$', fontsize = 28, y = 1.02)
	plt.plot((plt.xlim()[0], plt.xlim()[1]), (1, 1), 'k--')
	if config['fixedRange']:
		plt.ylim([0.8,2.])
	plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
	plt.tight_layout()
	plt.savefig(config['outputDir'] + 'METResponse_vs_BosonPt.pdf')
	plt.clf()
	plt.figure(10)
	plt.xlabel(r'$p_T^H\ \mathrm{in\ GeV}$', fontsize = 22)
	plt.ylabel(r'$\sigma(E_{T,||}^{miss}-E_{T,gen}^{miss})\ \mathrm{in\ GeV}$',fontsize = 20)
	plt.title(r'$\mathrm{Resolution}\ E_{T,||}^{miss}$', fontsize = 28, y = 1.02)
	legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
	#legend.get_frame().set_alpha(0.5)
	if config['fixedRange']:
		plt.ylim([16.,30.])
	plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
	plt.tight_layout()
	plt.savefig(config['outputDir'] + 'METResolution_vs_BosonPt.pdf')
	plt.clf()
	plt.figure(11)
	plt.xlabel(r'$p_T^H\ \mathrm{in\ GeV}$', fontsize = 22)
	plt.ylabel(r'$\sigma(E_{T,\bot}^{miss})\ \mathrm{in\ GeV}$',fontsize = 22)
	plt.title(r'$\mathrm{Resolution}\ E_{T,\bot}^{miss}$', fontsize = 28, y = 1.02)
	legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
	#legend.get_frame().set_alpha(0.5)
	if config['fixedRange']:
		plt.ylim([16.,30.])
	plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
	plt.tight_layout()
	plt.savefig(config['outputDir'] + 'METResolutionPerp_vs_BosonPt.pdf')
	plt.clf()
	plt.figure(12)
	legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
	#legend.get_frame().set_alpha(0.5)
	#plt.ylabel(r'$<U_{\|} / -p_T^H>$', fontsize = 20)
	plt.ylabel(r'$\left<-p_T^H/U_{||}\right>$', fontsize = 22)
	plt.xlabel(r'$p_T^H\ \mathrm{in\ GeV}$', fontsize = 22)
	plt.title(r'$\mathrm{Response^{-1}}\ U_{||}$', fontsize= 28, y = 1.02)
	plt.plot((plt.xlim()[0], plt.xlim()[1]), (1, 1), 'k--')
	if config['fixedRange']:
		plt.ylim([0.0,1.3])
	plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
	plt.tight_layout()
	plt.savefig(config['outputDir'] + 'ResponseInverted_vs_BosonPt.pdf')
	plt.clf()
	plt.figure(13)
	legend = plt.legend(loc='best', shadow=True, numpoints=1, fontsize = 20)
	#legend.get_frame().set_alpha(0.5)
	plt.xlabel(r'$p_T^H\ \mathrm{in\ GeV}$', fontsize = 22)
	plt.title(r'$\mathrm{Response^{-1}}\ E_{T,||}^{miss}$', fontsize= 28, y = 1.02)
	#plt.ylabel('MET_Long/genMET',fontsize = 20)
	plt.ylabel(r'$\left<E_{T,gen}^{miss}/E_{T,||}^{miss}\right>$',fontsize = 22)
	#plt.ylabel(r'$\ensuremath{{\not\mathrel{E}}_T}$',fontsize = 20)
	plt.plot((plt.xlim()[0], plt.xlim()[1]), (1, 1), 'k--')
	if config['fixedRange']:
		plt.ylim([0,2])
	plt.text(0.0*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0],1.015*(plt.ylim()[1]-plt.ylim()[0])+plt.ylim()[0],args.decayProcess,color = 'k',fontsize=18, weight='bold')
	plt.tight_layout()
	plt.savefig(config['outputDir'] + 'METResponseInverted_vs_BosonPt.pdf')
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
	else:
		print("Size of dataset: %i"%plotData.shape[0])

	print(plotData.shape)

	# perform plotting
	plot_results(config, plotData, dictPlot)

	#save the config as json in the plots folder
	with open(config['outputDir'] + 'config_%s.json'%os.path.basename(config['outputDir'][:-1]), 'w') as fp:
		json.dump(config, fp, sort_keys=True, indent=4)


	print('Plots created!')

	#exporting of the plots to ekpwww/nzaeh
	if 'export' in config:
		bashCommand = 'cp -r %s /ekpwww/web/nzaeh/public_html/'%config['outputDir']
		os.system(bashCommand)

		#copying index.php which allows viewing all plots at once on ekpwww by opening php.index
		bashCommand = 'cp index.php /ekpwww/web/nzaeh/public_html/%s/'%os.path.basename(config['outputDir'][:-1])
		os.system(bashCommand)

		print('Plots exported to ekpwww!')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Make MVAMet control plots.')
	parser.add_argument('-p', '--plottingconfig', default='../configs/config.json', help='Path to configurations file')
	parser.add_argument('-i', '--inputfile',nargs='+', default='', help='[optional] Inputfile(s) from which to create the plots from')
	parser.add_argument('-l', '--inputlabel',nargs='+', default='', help='[optional] Inputlabelname(s) to use in plots')
	parser.add_argument('-o', '--outputfolder', default='', help='[optional] Foldername in which to store the plots in')
	parser.add_argument('-d', '--decayProcess', default=r'$Z \to \mu\mu$', help='[optional] Decay process, Default: r"$Z \to \mu\mu$"')
	parser.add_argument('-c', '--constraints', default='', help='[optional] Constraints to data. E.g.: (50<=plotData[:,dictPlot["Boson_Pt"]]) & (50<=plotData[:,dictPlot["recoilslimmedMETs_LongZ"]]) NOTE: Use Single quotes in variable declaration')
	parser.add_argument('-e', '--export', dest='export', action='store_true', help='[optional] Exports plots to ekpwww after creating them')
	parser.add_argument('-s', '--study', dest='study', action='store_true', help='[optional] Create detailed parameter studies for each variable')
	parser.add_argument('-t', '--trueBoson', dest='trueBoson', action='store_true', help='[optional] Add gen Met to BosonPt to approximate BosonPt for genuine Met Events')
	parser.add_argument('-n', '--skip', dest='skipControlPlots', action='store_true', help='[optional] Skips creation of control plots, custom plots only')
	parser.add_argument('-r', '--fixedRange', dest='fixedRange', action='store_true', help='[optional] Fix ranges to default ranges in order to facilitate comparison between plots')
	parser.add_argument('-m', '--method',nargs='+', default=['Trunc'], help='[optional] Change method(s) to obtain Values. [Empiric, Trunc, Fit, FWHM, Quantile, Alpha, FWHMSpline, AlphaExcl, AlphaExclInt]')
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

	config['decayProcess'] = args.decayProcess
	print(config['decayProcess'])

	if args.fixedRange:
		config['fixedRange'] = True
		print('Ranges will be fixed to predefined ranges.')

	if args.study:
		config['studyVars'] = True
		print('Detailed variable studies will be performed.')

	if args.trueBoson:
		config['trueBoson'] = True
		print('Generator Met will be added to BosonPt to approximate true BosonPt for genuine Met events.')

	if args.skipControlPlots:
		config['skipControlPlots'] = True
		print('Only custom plots will be created.')

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


