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
import random
import math
import numpy as np

from . import chromModel
from . import trackModel
from . import SDDWriter

# Declare a bunch of constants
# # For proton, PDG code is 2212
# # For photon, PDG code is 22
DSBRate = 1.0
DSBComplexity = 0.43
directFrac = 0.4
DSBPerGy = 35

# Reference target radius for human cells, to calculate DNA density and derived energy per DSB
refTargRadius = 4.229
refEPerGy     = (4/3*math.pi*pow(refTargRadius,3))*6.242 # in keV, Using 1 Gy = 6.242 keV/um^3
refEPerDSB    = refEPerGy/DSBPerGy # For these values, this is 56.5 keV

# Should we use the sparse SDD format?
writeSparse=True

# Little handler to convert a list to comma separated string
def toCSV(A,sep=','):
	return sep.join(map(str,A))

# Return model damage string
def generateDmgandBase(damageType):
	fullBreakType = [0,0,0]

	damageArray = [[0*n for n in range(20)] for m in range(4)]
	# Single SSB
	if damageType==[0,1,0]:
		strand = random.choice([0,3])
		damageArray[strand][0]=1
		fullBreakType[1]+=1

	# Single DSB, losing 1 to 10 bases
	if damageType==[0,0,1]:
		breakLength = random.randint(1,10)
		# Single bp break is one column
		if breakLength==1:
			for row in range(4):
				damageArray[row][0]=1
		else:
			# Multi-BP break is marked by start and end
			for row in range(4):
				for position in range(breakLength):
					damageArray[row][0]=1
					damageArray[row][position+1]=1

		fullBreakType[0]=0
		fullBreakType[1]=2
		fullBreakType[2]=1

	# DSB+SSB within 10 BP
	if damageType==[0,1,1]:
		SSBGap = random.randint(-10,10)
		SSBStrand = random.choice([0,3])
		breakLength = random.randint(1,9)

		if SSBGap<0:
			SSB_BP = 0
			DSB_Start = SSBGap
		else:
			SSB_BP = SSBGap+breakLength
			DSB_Start = 0

		damageArray[SSBStrand][SSB_BP]=1
		fullBreakType[1]+=1

		if breakLength==1:
			for row in range(4):
				damageArray[row][DSB_Start]=1
		else:
			# Fill in all parts of a multi-BP entry
			for row in range(4):
				for position in range(breakLength):
					damageArray[row][DSB_Start+position]=1

		fullBreakType[0]=0
		fullBreakType[1]+=2
		fullBreakType[2]=1

	# Print out data. Row by row, left to right. 
	damageString = ''
	for row in range(4):
		for col in range(len(damageArray[row])):
			if damageArray[row][col]>0:
				damageString=damageString+' '.join(map(str,[row+1,col+1,damageArray[row][col]]))+'/'

	# Generate DNA bases string
	bases = [random.randint(1,4) for n in range(len(damageArray[0]))]

	# Remove single deletions
	for base in range(len(damageArray[0])):
		if damageArray[1][base]==1:
			bases[base]=0

	baseString = ''.join(map(str,bases))

	return damageString,baseString,fullBreakType

# Generate a damage uniformly distributed within cell
def XRayHits(DSBCount = 1.0,radius=1.0):
	# Poisson distribute around mean count number
	targetDSBs = np.random.poisson(DSBCount)
	# Also assign SSBs if calculating
	if DSBRate<1.0:
		SSBScaling = (1-DSBRate)/DSBRate
		targetSSBs = max(1,np.random.poisson(DSBCount/SSBScaling))
	else:
		targetSSBs = 0
	breakCounts = [targetSSBs,targetDSBs]

	# For each break, generate X,Y, Z positon randomly
	retBreaks = []
	for dsb,breaks in enumerate(breakCounts):
		for n in range(breaks):
			phi = 2*math.pi*random.random()
			theta = math.acos(random.uniform(-1,1))
			u = random.random()

			r = radius*pow(u,1.0/3.0)

			x = r*math.sin(theta)*math.cos(phi)
			y = r*math.sin(theta)*math.sin(phi)
			z = r*math.cos(theta)
			retBreaks.append( [x,y,z,dsb,1] )

	return retBreaks

# Generate damage distributed around an ion track
def ionHits(DSBCount=1.0, radius=1.0, LETdata=None, fixedTracks=None, breakStats=False):
	# If LET is 0, fall back to uniform X-ray distribution
	LET, radialData, energy, EScaling = LETdata
	if LET==0:
		return XRayHits(DSBCount,radius)

	# Calculate hits per track as it passes from -radius to +radius. 
	# Pad beam radius to target radius + 99% track radius
	DSBPerTrack = (LET/(refEPerDSB*EScaling))*(2.0*radius)
	padding = trackModel.sampleRadialPos(0.99,radialData)

	# Estimate mean number of tracks needed to deposit target number of hits in nucleus
	# Determined using ratio of nucleus volume to whole exposed cylinder
	trackEstimate = DSBCount/(DSBPerTrack*EScaling * (4/3*np.pi*radius*radius)/(2*np.pi*pow(radius+padding,2)))

	if fixedTracks is None:
		actualTracks = np.random.poisson(trackEstimate)
	else:
		actualTracks = fixedTracks

	# Generate damage by track
	retBreaks = []
	coreBreaks = 0
	rList = []

	for m in range(actualTracks):
		newEvent = 1

		# Calculate X,Y position where track arrives. Use target radius + 99% track radius
		u = random.random()
		r = (radius+padding)*pow(u,1.0/2.0)
		phi = 2*math.pi*random.random()
		X = r*math.cos(phi)
		Y = r*math.sin(phi)

		# Calculate actual DSB and SSB by this track
		trackDSBs = np.random.poisson(DSBPerTrack)
		if DSBRate<1.0:
			SSBScaling = (1-DSBRate)/DSBRate
			trackSSBs = np.random.poisson(DSBCount/SSBScaling)
		else:
			trackSSBs = 0
		breakCounts = [trackSSBs,trackDSBs]

		# For each break, position randomly along track length, and sample radial position from
		# track data file.
		for dsb,breaks in enumerate(breakCounts):
			for n in range(breaks):
				zPos = 2*(0.5-np.random.uniform())*radius

				radialFrac = np.random.uniform()
				dr = trackModel.sampleRadialPos(radialFrac,radialData)

				dPhi = 2*math.pi*random.random()
				xPos = X+dr*math.cos(dPhi)
				yPos = Y+dr*math.sin(dPhi)

				# Make sure that sampled hit actually remains within the nucleus
				if xPos*xPos + yPos*yPos + zPos*zPos < radius*radius:
					retBreaks.append([xPos,yPos,zPos,dsb,newEvent])
					newEvent = 0
					rList.append(dr)
					if dr<0.05:
						coreBreaks+=1

	# Print some statistics about break radial positions if requested
	if breakStats:
		print(LET, LET/(refEPerDSB*EScaling), actualTracks, len(retBreaks), coreBreaks, end=' ')
		print('\t'.join(map(str,(len([x for x in rList if x>r and x<r+0.0025])*1.0/len(rList) 
			                          for r in np.arange(0,4,0.0025)))))
	return retBreaks

# Take a list of positions and complexities, and format them into full breaks
def formatBreaks(breakPositions,radius=1.0, bdRange=-1, letData=None, particleTypes="2212", 
	             timeProfile = None, firstField = True):
	newHits = []
	eventNo=0
	# For each position, generate break characteristics
	for pos in breakPositions:
		x,y,z,dsb,newEvent = pos

		# Increment counter if it's a new event and set particle time
		if newEvent>0: 
			eventNo+=1
			if timeProfile==None:
				pTime = 0
			else:
				pTime = timeProfile[0] + random.random()*timeProfile[1]

		# If hit list is empty and this is the first field, this is a new exposure
		if firstField and newHits==[]:
			newEvent = 2
			eventNo = 0 #Reset to zero for new exposure

		# Build 3D break extent
		breakExtent  = [x,y,z,x+0.01,y+0.01,z+0.01,x-0.01,y-0.01,z-0.01]
		extentString = (toCSV(breakExtent[0:3])+ '/'+toCSV(breakExtent[3:6])+'/'+
					    toCSV(breakExtent[6:9]) )

		# Set break type
		if dsb==0:
			breakType = [0,1,0]
		else:
			if random.random()>DSBComplexity:
				breakType = [0,0,1]
			else:
				breakType = [0,1,1]

		# Sample chromosome, and generate illustrative damage structure
		chromID, chromPos = chromModel.modelChromosome(x,y,z)
		damageString,baseString,fullBreakType = generateDmgandBase(breakType)

		# Set BDs to 0 if we're not logging those
		if bdRange<0: fullBreakType[0]=0

		# Set cause through random assignment - direct or indirect, all breaks
		if random.random()>directFrac:
			cause = [1, 0, sum(fullBreakType[:-1])]
		else:
			cause = [0, sum(fullBreakType[:-1]), 0]

		# Placeholder values for other parameters if a full output is requested
		time = str(random.random()*2.0)
		energies = letData[2]
		trans = [x,y,-radius]
		direction = [0,0,0]

		# Append hit data, either in minimal or comprehensive format
		if writeSparse:
			newHits.append([toCSV([newEvent,eventNo],','), toCSV([x,y,z],', '), 
								toCSV(fullBreakType)])
		else:
			newHits.append([toCSV([newEvent,eventNo],','),extentString, chromID, 
				                chromPos, toCSV(cause), toCSV(fullBreakType), damageString, baseString,
				                time, particleTypes, energies, toCSV(trans,'/'), 
				                toCSV(direction,'/'), pTime])
	return newHits

def simFromFile(posFile, chromosomes=46, letData = [1, None, 1.0], incident ="2212", dose=-1):
	with open(posFile) as inFile:
		# Get initial radius and prep average chromosomes
		inputRow = inFile.readline()
		while inputRow[0]=='#': inputRow = inFile.readline()
		radiusData = [float(r) for r in inputRow.strip().split('\t')]
		meanRadius = pow(np.product(radiusData),1/len(radiusData))
		chromModel.subDivideSphere(chromosomes,meanRadius)

		allBreaks = []
		breakPositions = []
		for row in inFile:
			row = row.strip()
			# Blank lines delimit new exposures
			if len(row)==0:
				if len(breakPositions)>0: allBreaks.append(breakPositions)
				breakPositions = []
				continue

			# Get break position
			breakPos = [float(r) for r in row.split('\t')]
			# Append a 1 to denote a DSB, and flag new event as appropriate
			breakPos.append(1)
			if len(breakPositions)==0:
				breakPos.append(2)
			else:
				breakPos.append(0)
			breakPositions.append(breakPos)
		if len(breakPositions)>0: allBreaks.append(breakPositions)

	outBreaks = []
	for breakPositions in allBreaks:
		outBreaks.append(formatBreaks(breakPositions, letData = letData))

	outName = posFile.split('.')[0]+'.sdd'
	SDDWriter.writeToFile(outBreaks,outName, writeSparse, meanRadius, geometry = [1]+radiusData, 
						  incident=incident, dose=dose)


# Generate hits for a given exposure
def generateHits(runs=1, radius=1.0, DSBCount=1, chromosomes=1, bdRange=-1, letData=None, 
				 particleTypes="2212", timeProfile = None, firstField = True):
	# If this is an entirely new exposure, update chromosome model
	if firstField: chromModel.subDivideSphere(chromosomes,radius)
	hitList = []

	for n in range(runs):
		breakPositions = ionHits(DSBCount, radius, letData)
		hitList.append(formatBreaks(breakPositions, radius, bdRange, letData, particleTypes, 
									timeProfile, firstField))

	return hitList

# Generate hits for a requested exposure, and write this data to an SDD-formatted file
def simExposure(hits, runs, chromosomes, outFile, targetVol, geometry=[1,3,3,3], DNADensity=-1,
	            bdRange=-1, O2=-1, incident="22", energy=0.1, function="Point", 
	            grouping="Single Event", letData=None, timeProfile = None):
	hitData = generateHits(runs, geometry[1], hits, chromosomes, bdRange, 
						   letData, particleTypes=incident, timeProfile = timeProfile)
	SDDWriter.writeToFile(hitData, outFile, writeSparse, targetVol, geometry, DNADensity, bdRange,
						  O2, incident, energy, hits/DSBPerGy, function, grouping)

# Generate hits for a multi-field exposure, and write this data to an SDD-formatted file.
# Syntax same as single exposure above, but relevant parameters: hits, incident, letData, and
# timeProfile should be passed in as lists
def simMultiExposure(hits, runs, chromosomes, outFile, targetVol, geometry=[1,3,3,3], DNADensity=-1,
	              bdRange=-1, O2=-1, incident="22", energies=0.1, function="Point", 
	              grouping="Single Event", letData=None, timeProfile = None):
	hitData = None
	for h,let, particle, times in zip(hits,letData,incident,timeProfile):
		# Just generate hits for first field
		if hitData == None:
			hitData = generateHits(runs, geometry[1], h, chromosomes, bdRange, let, 
								   particleTypes=particle, timeProfile = times)
		else:
		# For later fields, need to merge into each run
			newHits = generateHits(runs, geometry[1], h, chromosomes, bdRange, let, 
								   particleTypes=particle, timeProfile = times, firstField = False)
			for n in range(len(hitData)):
				hitData[n]+=newHits[n]

	# Make sure first event in each exposure, assuming something occurs is marked as a new exposure (2)
	for n in range(len(hitData)): 
		if len(hitData[n])>0:
			hitData[n][0][0]='2'+hitData[n][0][0][1:]

	doses = toCSV([h/DSBPerGy for h in hits])
	SDDWriter.writeToFile(hitData, outFile, writeSparse, targetVol, geometry, DNADensity, bdRange,
						  O2, incident, energies, doses, function, grouping)

# Generate PID for different particle types, based on atomic number (Z=0 is gamma)
def PIDLookup(Z):
	if Z==0:
		return '22'
	if Z==1:
		return '2212'
	if Z>1:
		# Base string - nucleus with 0 strange quarks
		baseString = '100'
		zVal = str(Z)
		while len(zVal)<3: zVal = '0'+zVal
		aVal = str(2*Z)
		while len(aVal)<3: aVal = '0'+aVal

		return baseString+zVal+aVal+'0'

# Look up appropriate datafile, if available. Fall back to Carbon if non-supported ion is requested
def dataFileNames(Z):
	particleNames = ['Gamma', 'Proton', 'Helium', 'Lithium', 'Beryllium', 'Boron', 'Carbon', 
					 'Nitrogen', 'Oxygen']
	if Z==1 or Z==2 or Z==6 or Z==7:
		return 'Radial Energy '+particleNames[Z]+'.xlsx'
	else:
		return 'Radial Energy Carbon.xlsx'

# Generate an exposure for a requested ion, LET, and number of repeats
def generateExposure(energy, LET, dose, particleZ, runs, targetRadius=4.32, 
				     chromosomes=46, timeProfile = None, extraTargetInfo='',
				     fileName = None):
	# Some general model parameters for SDD header
	DNADensity = 6100/(4.0/3.0*math.pi*pow(refTargRadius,3))
	geometry = [1,targetRadius,targetRadius,targetRadius]
	targetVol = str(targetRadius)+' um spherical nucleus'
	particleID = PIDLookup(particleZ)
	bdRange = -1 # Don't record base damages

	# Set up output file name
	# Calculate hits by scaling by dose
	hits = dose*DSBPerGy

	if fileName ==None:
		fileBase = 'DNA Damage Z='+str(particleZ)+ ' '
		fileName = fileBase+str(energy)+' MeV '+ str(dose) + ' Gy'
		if writeSparse:
			fileName = fileName + ' minimal.txt'
		else:
			fileName = fileName + ' full.txt'

	# Get track data for model input
	letData = None
	if particleZ>0:
		trackFile = dataFileNames(particleZ)
		trackModel.readCumuDoseFile(trackFile)
		letData = [LET, trackModel.buildCumCurve(LET), energy, pow(targetRadius/refTargRadius,3)]
	else:
		letData = [0,None,energy,pow(targetRadius/refTargRadius,3)]

	# Run simulation
	print(energy, LET, hits)
	simExposure(hits, runs, chromosomes, fileName, targetVol, geometry, DNADensity, bdRange,
				letData=letData, incident=particleID, energy=energy, timeProfile=timeProfile)

# Generate an exposure for an abtirary set of exposures
def generateMultiExposure(energies, LETs, doses, particleZs, runs, targetRadius=4.32, 
				          chromosomes=46, timeProfiles = None, extraTargetInfo='', 
				          fileName = None):
	# Some general model parameters for SDD header
	DNADensity = 6100/(4.0/3.0*math.pi*pow(refTargRadius,3))
	geometry = [1,targetRadius,targetRadius,targetRadius]
	targetVol = str(targetRadius)+' um spherical nucleus'
	particleID = [PIDLookup(z) for z in particleZs]
	bdRange = -1 # Don't record base damages

	# Set up output file name
	# Calculate hits by scaling by dose
	hits = [d*DSBPerGy for d in doses]

	if fileName ==None:
		fileBase = 'DNA Damage Zs '+toCSV(particleZs, ' ')+' '
		fileName = fileBase+toCSV(energies, ' ')+' MeV '+ toCSV(doses, ' ') + ' Gy'
		if writeSparse:
			fileName = fileName + ' minimal.txt'
		else:
			fileName = fileName + ' full.txt'

	# Get track data for model input
	letData = []
	for particleZ, LET, energy in zip(particleZs, LETs, energies):
		if particleZ>0:
			trackFile = dataFileNames(particleZ)
			trackModel.readCumuDoseFile(trackFile)
			letData.append([LET, trackModel.buildCumCurve(LET), energy, 
							pow(targetRadius/refTargRadius,3)] )
		else:
			letData.append([0, None, energy,
							pow(targetRadius/refTargRadius,3)] )

	# Run simulation
	print(energies, LETs, hits)
	print(fileName)
	simMultiExposure(hits, runs, chromosomes, fileName, targetVol, geometry, DNADensity, bdRange,
				     letData=letData, incident=particleID, energies=energies, timeProfile=timeProfiles)


# Example application
# Build a basic X-ray and ion dataset
def basicXandIon(targetRadius = 4.32, runs = 10, conditions=None, extraTargetInfo = ''):
	# Photons, doses from 1 to 8 Gy
	particleZ = 0
	ionConditions = [[1.0,0,1], [1.0,0,2], [1.0,0,3], [1.0,0,4], 
					 [1.0,0,6], [1.0,0,8]]
	for energy, LET, dose in ionConditions:
		generateExposure(energy, LET, dose, particleZ, runs, targetRadius, chromosomes=46,
						 timeProfile = [0,60*1E9*2])

	# Protons, at a dose of 1 Gy 
	particleZ = 1 
	ionConditions = [ [0.975,29.78,1], [1.175,25.27,1], [1.5,20.59,1], 
				      [1.8, 17.78, 1], [2.2,15.19,1],   [2.5,13.72,1], 
				      [3.5, 10.60, 1], [5.5,7.42,1],    [8.5,5.25,1], 
				      [34,  1.77,  1] ]
	for energy, LET, dose in ionConditions:
		generateExposure(energy, LET, dose, particleZ, runs, targetRadius, chromosomes=46,
						 timeProfile = [0,60*1E9*2])

	# Carbon ions, at a dose of 1 Gy
	particleZ = 6
	ionConditions = [ [24,512,    1], [60,265,1], [120,151.95,1],
					  [185,100,   1], [360,60,1], [960,26,1],
					  [1200,20.29,1] ]
	for energy, LET, dose in ionConditions:
		generateExposure(energy, LET, dose, particleZ, runs, targetRadius, chromosomes=46,
						 timeProfile = [0,60*1E9*2])
		