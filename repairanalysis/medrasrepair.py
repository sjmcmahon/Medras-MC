# ############################################################################
# 
# This software is made freely available in accordance with the simplifed BSD
# license:
# 
# Copyright (c) <2018>, <Stephen McMahon>
# All rights reserved
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation 
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
# DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contacts: Stephen McMahon,	stephen.mcmahon@qub.ac.uk
# 
# ############################################################################
#
# Generic Medras-MC wrapper. Provides full fidelity output, misrepair separation, and misrepair 
# spectrum output through different functions.
#
# ############################################################################
import numpy as np
import os
import re
import copy

from . import medrasparser
from . import misrepaircalculator as calcMR 
from . import analyzeAberrations

# Input parameters common to different processes
sigma = 0.04187 # Misrepair range, as fraction of nuclear radius
maxExposures = 1000  # Maximum exposures per file to simulate
repeats = 50		 # Number of repeats for each exposure
minMisrepSize = 0 	 # Neglect misrepair events separated by less than this genetic distance (MBP)

# Options for fidelity repair output
writeKinetics = True
writeAllKinetics = False
addFociDelay = True
kineticLimit = 25 	 # Hours, maximum time for repair kinetics

# Options for misrepair spectrum output
doPlot = False
allFragments = False
listAcentrics = False
simulationLimit = np.inf # Hours, time at which to simulate misrepair

# Method to sort files nicely with numbers.
def sort_nicely( l ):
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
	l.sort( key=alphanum_key )

###################
#
# Misrepair spectrum analysis and helpers
# 
###################
def prepareDamage(misrepairList, remainingBreaks, chromosomes):
	# Trim misrepair list - store chromosome, genetic position, and orientation of each break (a&b)
	trimMisrep = [ [[a[3][1],a[4],a[5],a[1]], [b[3][1],b[4],b[5],a[1]] ] for a,b,c,d in misrepairList] 
	# Extract same data for remaining breaks
	trimBreaks = [ [a[3][1],a[4],a[5],a[1]] for a in remainingBreaks]

	# Scale by chromosome sizes, if needed
	for misrep in trimMisrep:
		for damage in misrep:
			chromID = damage[0]
			if damage[1]<1: damage[1] *= chromosomes[chromID][2]
	for damage in trimBreaks:
		chromID = damage[0]
		if damage[1]<1: damage[1] *= chromosomes[chromID][2]

	# Filter out duplicated DNA ends
	for n in range(len(trimBreaks)-1,0,-1):
		if trimBreaks[n][0]==trimBreaks[n-1][0] and trimBreaks[n][1]==trimBreaks[n-1][1]:
			trimBreaks.pop(n)

	return trimMisrep, trimBreaks

def listAcentricSizes(baseChromosomes, finalChromosomes,pos=0.5):
	fullChroms = analyzeAberrations.characteriseChroms(finalChromosomes)
	print('\n\nAcentric sizes:')
	for c in fullChroms:
		centromereCount = analyzeAberrations.centricCount(c,baseChromosomes,pos)

		if centromereCount==0:
			theSize = sum( abs(f[1]-f[2]) for f in c[3])
			print(theSize)

def misrepairSpectrum(fileData, header, fileName):
	allBreaks, scaledSigma, meanE, complexity, complexitySd, dose, emptySets = fileData
	radius = scaledSigma/sigma

	print('Data for: '+fileName)
	if len(allBreaks)==0:
		print('No damage in file')
		return None

	# Test to see if chromosome ID data is available, abort if not
	firstBreak = allBreaks[0][0]
	if firstBreak[3][0]==-1:
		print('No chromosome IDs found in data file, aborting!')
		return None

	noChroms = header['Chromosomes'][0]
	baseChromosomes = [[n,0,header['Chromosomes'][1][n] ] for n in range(noChroms) ]
	fullFrags = []
	allChroms = []
	allRings = []
	for m,breakList in enumerate(allBreaks):
		if m>=maxExposures:
			break
		misrepList,repList, remBreaks = calcMR.singleRepair(copy.deepcopy(breakList), None, 
															scaledSigma, finalTime = simulationLimit)
		trimMisrep, trimRemBreaks = prepareDamage(misrepList, remBreaks, baseChromosomes)

		chroms, rings, frags = analyzeAberrations.doRepair(baseChromosomes, trimMisrep, 
			                        remBreaks = trimRemBreaks, index=m, breaks=len(breakList)//2, 
			                        baseBreaks=breakList, plot = doPlot, allFragments=allFragments, 
			                        outFile = fileName+str(m)+'.png')

		frags = [f+[m] for f in frags]
		fullFrags += frags
		allChroms = allChroms + chroms 
		allRings = allRings + rings 
	if listAcentrics:
		listAcentricSizes(baseChromosomes, allChroms+allRings)

	return None

###################
#
# Fidelity analysis and helpers
# 
###################

# Convert a list of repair times to a normalised kinetic curve
def summariseKinetics(repairTimes):
	repairTimes.sort()
	timeStep = 0.1
	currTime = 0
	kinetic = [1]
	for n, t in enumerate(repairTimes):
		while t>currTime+timeStep:
			kinetic.append(1.0-1.0*n/len(repairTimes))
			currTime+=timeStep
			if currTime+timeStep>kineticLimit:
				break
		if currTime+timeStep>kineticLimit:
			break
	kinetic = kinetic+[0]*(int(kineticLimit/timeStep)-len(kinetic))
	return '\t'.join(map(str,kinetic))

# Generate a set of summary statistics for a file
def summariseFidelity(fileName, complexity, outputs):
	totalBreaks = 1.0*sum(o[0] for o in outputs)
	averageBreaks = np.mean([o[0] for o in outputs])
	breakStdev = np.std([o[0] for o in outputs])

	if totalBreaks>0:
		averageMisrep = np.mean([o[0]*o[1] for o in outputs])/averageBreaks
		misrepStdev   = np.std([o[1] for o in outputs])
		fileAverages = [sum([o[n]*o[0] for o in outputs])/totalBreaks for n in range(3,6)]
	else:
		averageMisrep=0 
		misrepStdev = 0 
		fileAverages = [0,0,0]
	smry = (fileName+'\tSummary\t'+str(totalBreaks)+'\t'+str(complexity)+'\t'+str(averageBreaks)+
		   '\t'+str(breakStdev)+'\t'+str(averageMisrep)+'\t'+str(misrepStdev)+'\t'+
		   str(fileAverages[0]) )
	#smry+= '\t\t'+'\t'.join(map(str,fileAverages[1:]))
	return smry

fidelityRun = False
def repairFidelity(fileData, header, fileName):
	# Print output header
	print('File\tBreak Set\tBreak Count\tMisrepair\tStdev\tInter-Chromosome Rate', end='')
	#print('\t\tAnalytic Misrepair\tEtaSum', end='')
	if writeAllKinetics:
		print('\t\t', '\t'.join(map(str,[tau/10.0 for tau in range(10*kineticLimit)])), end='')
	print()

	# Extract file data
	breaks, scaledSigma, meanE, complexity, complexitySd, dose, emptySets = fileData
	radius = scaledSigma/sigma
	chromSizes = header['Chromosomes'][1]

	outputs = [[0,0,0,0,0,0,0]]*emptySets
	outputTimes = []
	# Analyze each set of breaks
	for n,breakList in enumerate(breaks):
		if n>=maxExposures: break
		print(fileName+'\t'+str(n),'\t', len(breakList)/2,'\t', end='')

		if len(breakList)<=2 or len(breakList)>20000:
			repairData = [0 ,0 ,0, [], 0, 0]
			if len(breakList)>20000:
				print ('Skipping due to memory concerns','\t', end= '')
		else:
			repairData = calcMR.fullRepair(breakList, scaledSigma, repeats=repeats, 
										   addFociClearance=addFociDelay, radius=radius, 
										   chromSizes = chromSizes, sizeLimit = minMisrepSize)

		# Unpack repair data
		mean, stdevRate, interChromRate, repTimes, analyticMisrep, etaSum = repairData
		outputs.append([len(breakList)/2.0, mean, stdevRate, 
					   interChromRate, analyticMisrep, etaSum])
		outputTimes+=repTimes

		# Write misrepair and kinetic data
		print(mean,'\t', stdevRate,'\t',interChromRate, end='')
		#print('\t\t','\t'.join(map(str,[analyticMisrep, etaSum])), end='')
		if writeAllKinetics:
			print('\t\t', summariseKinetics(repTimes), end='')
		print()
	print()

	# Build summary of file
	summary = ''
	global fidelityRun
	if fidelityRun is False:
		fidelityRun = True
		summary += ('\nFile\tBreak Set\tTotal Breaks\tComplexity\tBreaks per Exposure\tBreaks Stdev'
			        '\tMisrepair\tStdev\tInter-Chromosome Rate')
		# summary+='\t\tAnalytic Misrepair\tEtaSum'
		if writeKinetics:
			summary += ('\t\tTime\t'+'\t'.join(str(tau/10.0) for tau in range(10*kineticLimit)))
		summary += '\n'

	summary+=summariseFidelity(fileName,complexity,outputs)
	if writeKinetics:
		summary+='\t\t'+fileName+'\t'
		summary+=summariseKinetics(outputTimes)

	return summary

###################
#
# Misrepair separation calculator
# 
###################
separationRun = False
def misrepairSeparation(fileData,header,fileName):
	# Initialise some shared values
	bins = 500
	breaks, scaledSigma, meanE, complexity, complexitySd, dose, emptySets = fileData
	radius = scaledSigma/sigma
	maxSeparation = 2*radius
	rBins = [n*(maxSeparation*1.0/bins) for n in range(bins+1)]

	# For first run, print header
	global separationRun
	if separationRun is False:
		separationRun = True
		print('File\tBreakSets\tEmptySets\tBreaks per set\tMisrepairs per set\tSeparation (um):\t', end=' ')
		print('\t'.join(map(str,rBins)))

	print(fileName,'\t',len(breaks),'\t',emptySets,'\t', end=' ')
	misrepairSeps = []
	totBreaks = 0
	m=0
	for m,breakList in enumerate(breaks):
		if m>=maxExposures:
			break
		totBreaks+=len(breakList)
		for n in range(repeats):
			misrepList,repList, remBreaks = calcMR.singleRepair(copy.deepcopy(breakList), None, scaledSigma)
			misrepairSeps+=[mis[2] for mis in misrepList]

	m = m+1 # Total exposures counted
	if len(misrepairSeps)>0:
		print(totBreaks/m/2, '\t',len(misrepairSeps)/(m*repeats),end='\t')
		print(fileName,'\t','\t'.join(map(str, np.histogram(misrepairSeps, rBins, density=True)[0]))) 
	else:
		print()
	return None

###################
#
# DSB separation calculator
# 
###################
def dsbSeparation(fileData,header,fileName):
	# Initialise some shared values
	bins = 500
	breaks, scaledSigma, meanE, complexity, complexitySd, dose, emptySets = fileData
	radius = scaledSigma/sigma
	maxSeparation = 2*radius
	rBins = [n*(maxSeparation*1.0/bins) for n in range(bins+1)]

	# For first run, print header
	global separationRun
	if separationRun is False:
		separationRun = True
		print('File\tBreakSets\tEmptySets\tSeparation (um):\t', end=' ')
		print('\t'.join(map(str,rBins)))

	print(fileName,'\t',len(breaks),'\t',emptySets,'\t',fileName,'\t', end=' ')
	seps = []
	for m, breakList in enumerate(breaks):
		#print(m,breakList)
		for i in range(0,len(breakList),2):
			for j in range(i+2,len(breakList),2):
				seps.append(np.sqrt(calcMR.distanceToSq(breakList[i][1],breakList[j][1])))
				#if np.sqrt(calcMR.distanceToSq(breakList[i][1],breakList[j][1]))<0.01:
				#	print(i,j)

	if len(seps)>0:
		print('\t'.join(map(str, np.histogram(seps, rBins, density=True)[0]))) 
	else:
		print()	

	return None

###################
#
# Radial damage separation calculator
# 
###################
radialRun = False
def radialDSBs(fileData,header,fileName):
	# Initialise some shared values
	bins = 1600
	breaks, scaledSigma, meanE, complexity, complexitySd, dose, emptySets = fileData
	radius = scaledSigma/sigma
	maxSeparation = 4.0
	rBins = [n*(maxSeparation*1.0/bins) for n in range(bins+1)]

	# For first run, print header
	global radialRun
	if radialRun is False:
		radialRun = True
		print('File\tBreakSets\tEmptySets\tSeparation (um):\t', end=' ')
		print('\t'.join(map(str,rBins)))

	print(fileName,'\t',len(breaks),'\t',emptySets,'\t',fileName,'\t', end=' ')
	seps = []
	for m, breakList in enumerate(breaks):
		for i in range(len(breakList)):
			#print('\t'.join(map(str,breakList[i][1])))
			seps.append(np.sqrt(pow(breakList[i][1][0],2)+pow(breakList[i][1][1],2)))
	if len(seps)>0:
		print('\t'.join(map(str, np.histogram(seps, rBins, density=True)[0]))) 
	else:
		print()	

	return None

###################
#
# Output damage by track
# 
###################
trackRun = False
def trackBreaks(fileData,header,fileName):
	# Initialise some shared values
	breaks, scaledSigma, meanE, complexity, complexitySd, dose, emptySets = fileData
	radius = scaledSigma/sigma

	# For first run, print header
	global radialRun
	if radialRun is False:
		radialRun = True
		print('File\tBreak count\tBreaks')

	totalTracks = 0
	totalBreaks = 0
	for m, breakList in enumerate(breaks):
		breakCount = 0
		for i in range(len(breakList)):
			if breakList[i][0]!=lastTrack:
				if lastTrack!=-1:
					print(fileName,lastTrack,breakCount)
					totalTracks+=1
					totalBreaks+=breakCount
				lastTrack = breakList[i][0]
				breakCount = 1
			else:
				breakCount+=1

		if lastTrack!=-1:
			print(fileName,lastTrack,breakCount)

	retString = "\t".join(map(str,[fileName,totalTracks,totalBreaks,totalBreaks/totalTracks]))
	return retString

###################
#
# Simulation wrapper
# 
###################
def repairSimulation(folder, analysisFunction='Fidelity', verbose=False):
	functions = [ ['Fidelity',  repairFidelity],
				  ['Separation',misrepairSeparation],
				  ['Spectrum',  misrepairSpectrum],
				  ['DSBSeparation',dsbSeparation],
				  ['DSBRadial',    radialDSBs],
				  ['TrackBreaks',  trackBreaks] ]
	summaries=[]
	timeSummary = []

	if folder[-1]!='/': folder=folder+'/'
	print('From folder:\t',folder,'\n')
	fileNames = os.listdir(folder)
	sort_nicely(fileNames)
	filePaths = [folder+f for f in fileNames]
	for filePath in filePaths:
		if os.path.isdir(filePath) or (filePath[-3:]!='txt' and filePath[-3:]!='sdd'):
			continue

		# Parse file
		fileName = filePath.split('/')[-1]
		fileData = medrasparser.parseToBreaks(filePath,sigma, verbose=verbose)
		header = fileData[-1]

		# Find matching method and run
		for name, func in functions:
			if analysisFunction == name:
				summaries.append(func(fileData[:-1], header, fileName))

		if len(summaries)==0:
			print('Did not find matching analysis function. Options are:')
			print(', '.join(f[0] for f in functions))
			return

	if len(summaries)==0:
		print('No output returned!')
		return

	if summaries[0] is not None:
		for summary in summaries:
			print(summary)
