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

import copy
import numpy as np

from . import plotAberrations

minFragment = 0   # MBP, minimum measurable fragment for DNA loss
maxFragment = 5.7 # MBP, maximum measurable fragment for DNA loss
largeMisrepThreshold = 3.0 # MBP, Size of 'large' misrepair

#######################################################
#
# Methods to calculate statistics on groups of misrepair
#
#######################################################
# Total misrepair, large misrepair, and inter-chromosome misrepair
def misrepairStats(misreps, chromosomes):
	interChrom = 0
	largeMisrep = 0
	for misrep in misreps:
		if misrep[0][0]!=misrep[1][0]:
			interChrom+=1
		else:
			chromID = misrep[0][0]
			chromLength = chromosomes[chromID][2]
			if abs((chromLength*(misrep[0][1]-misrep[1][1])) )>=largeMisrepThreshold:
				largeMisrep+=1
	return str(len(misreps))+"\t"+str(largeMisrep)+"\t"+str(interChrom)

# Fragment complexity - complex breaks are those with >1 chromosome junction
def calculateComplexities(finalFragments):
	simple = 0
	complexes = 0
	for f in finalFragments:
		junctions = 0
		for n in range(1,len(f)):
			if(f[n-1][0]//2!=f[n][0]//2):
				junctions+=1

		if junctions==1:
			simple+=1
		else:
			if junctions>1:
				complexes+=1

	return str(simple)+"\t"+str(complexes)

# Count centromeres in a single fragment
# Centromeres are defined at a fractional position given by pos
def centricCount(c, baseChromosomes, pos=0.5):
	centromereCount=0
	for f in c[3]:
		start = min(f[1],f[2])
		end = max(f[1],f[2])
		cent = pos*baseChromosomes[f[0]][2]
		if cent>=start and cent<=end:
			centromereCount+=1
	return centromereCount

# Yield of normal, acentric and multi-centric fragments
def centricCheck(baseChromosomes, finalChromosomes,pos=0.5):
	normal = 0
	acentric = 0
	multicentric = 0
	largeLoss = 0
	for c in finalChromosomes:
		centromereCount = centricCount(c,baseChromosomes,pos)

		if centromereCount==0:
			acentric+=1
			if sum( abs(f[1]-f[2]) for f in c[3])>3:
				largeLoss+=1
		if centromereCount==1:
			normal+=1
		if centromereCount>1:
			multicentric+=1
	return '\t'.join(map(str,[normal,acentric,multicentric,largeLoss]))

# Calculate initial fragment distribution
def fragmentDistribution(baseChromosomes,breakList):
	breakPoints = []
	for damage in breakList[::2]:
		breakPoints.append([damage[3][1],damage[4]*baseChromosomes[damage[3][1]][2]])

	# Sort by chromosome and genetic position
	breakPoints.sort(key=lambda x:(x[0],x[1]))
	currChrom = breakPoints[0][0]
	currPos = 0
	totalFrags = 0

	for damage in breakPoints:
		if currChrom!=damage[0]:
			fullChromSize = baseChromosomes[currChrom][2]
			remLength = fullChromSize-currPos
			if remLength<maxFragment and remLength>minFragment:
				totalFrags+=remLength
			currPos = 0
			currChrom=damage[0]
		fragLength = damage[1]-currPos
		if fragLength<maxFragment and fragLength>minFragment:
			totalFrags+=fragLength
		currPos = damage[1]
	return totalFrags

# Estimate of DNA lost in acentric fragments
def dnaLoss(baseChromosomes,finalChromosomes):
	lostDNA = 0
	allFragments = []
	for c in finalChromosomes:
		length = c[0]
		spatialSep = c[-1]
		if length<maxFragment:
			lostDNA+=length
		allFragments.append([length,spatialSep,c[1]])
	return lostDNA, allFragments

# Identify size of each chromosome, major component and if it's a ring
def characteriseChroms(chroms,rings=False,doPrint=False):
	retChroms = []
	for c in chroms:
		totLen = 0
		lenDict = dict()
		for f in c:
			totLen +=abs(f[2]-f[1])
			if f[0] not in lenDict:
				lenDict[f[0]]=abs(f[2]-f[1])
			else:
				lenDict[f[0]]+=abs(f[2]-f[1])

			if doPrint:
				print(f[0],abs(f[2]-f[1]),'\t', end=' ')

		if c[0][3] is None or c[-1][4] is None:
			spatialLen = -1
		else:
			spatialLen = np.linalg.norm(np.array(c[0][3])-np.array(c[-1][4]))

		if doPrint:
			print('\t',totLen)
		mainChrom = max(iter(lenDict.keys()), key=(lambda key: lenDict[key]))
		retChroms.append([totLen,mainChrom,rings,c,spatialLen])
	return retChroms


#######################################################
#
# Methods to build list of rejoined chromosomes
#
#######################################################

# Split initial chromosomes at sites of breaks
def splitChromosomes(chromosomes,breaks):
	chromStack = copy.deepcopy(chromosomes)
	chromStack = [ [c+[None,None]] for c in chromStack]
	for b in breaks:
		breakChrom = b[0]
		breakPos = b[1]
		spatialPos = b[3]
		for n,chrom in enumerate(chromStack):
			if chrom[0][0]==breakChrom and breakPos>chrom[0][1] and breakPos<chrom[0][2]:
				chromStack.pop(n)
				topChrom = [chrom[0][0],chrom[0][1],breakPos,chrom[0][3],spatialPos]
				chromStack.append([topChrom])
				botChrom = [chrom[0][0],breakPos,chrom[0][2],spatialPos,chrom[0][4]]
				chromStack.append([botChrom])
	return chromStack

# Locate chromosome in list matching rejoining event
def indexChrom(chromFragments, c, p, d):
	# Tolerate match within 1 part in 10^12, in case of rounding issues
	for n, chrom in enumerate(chromFragments):
		firstChrom = chrom[0]
		if firstChrom[0]==c:
			if abs(firstChrom[1]-p)<1E-12*p:
				if ( (firstChrom[2]-firstChrom[1]>0 and d>0) or 
				    (firstChrom[2]-firstChrom[1]<0 and d<0) ):
					return n,0
		lastChrom = chrom[-1]
		if lastChrom[0]==c:
			if abs(lastChrom[2]-p)<1E-12*p:
				if ( (lastChrom[2]-lastChrom[1]>0 and d<0) or 
				     (lastChrom[2]-lastChrom[1]<0 and d>0) ):
					return n,1
	print('Failed to find matching chromosome fragment! Requested chromosome details:')
	print(c,p,d)
	print('All Fragments:')
	for c in chromFragments:
		print(c)
	return None

# Append chromosome fragments in appropriate order, potentially inverting
def appendFragments(chromOne,chromTwo,end1,end2):
	# Need to flip one chromosome if it's RHS-RHS or LHS-LHS joining
	if end1==end2:
		chromTwo = [ [c[0],c[2],c[1],c[4],c[3]] for c in reversed(chromTwo)]
	# Line up appropriately
	if end1==0:
		return chromTwo+chromOne
	else:
		return chromOne+chromTwo

# Check if header needs to be printed on first execution
def checkHeader(initialBreaks=None, headerPrint=[False]):
	# Use mutable arguments to check if run before
	if headerPrint[0]:
		return
	headerPrint[0] = True
	print('Index\tBreaks\tResidual\tMisrepairs\tLarge Misrepairs\tInter-Chromosome Misrepairs\t'
		  'Single-Junction Chromosomes\tMulti-Junction Chromosomes\tNormal Chromosomes\t'
		  'Acentric Linear\tMulti-Centric\tLarge Loss\tCentric Ring\tAcentric Ring\t'
		  'Multi-Centric Ring\tLarge Loss', end='')
	if initialBreaks is not None:
		print('\tInitial DNA Fragmentation\tPotential DNA loss', end='')
	print()

# Core repair loop, iterate over each repair and append chromosome fragments
def doRepair(chromosomes, repairs, remBreaks=None, index=0, breaks=-1, baseBreaks=None, plot=False, 
			 allFragments=False, inFile=None, outFile=None):
	checkHeader(baseBreaks)
	# Build our breaklist, and abort if empty
	breakList = [b for r in repairs for b in r]
	if remBreaks is not None:
		breakList = breakList+remBreaks
	if len(breakList)==0:
		print(index,'\tNo misrepair!')
		return [],[], []

	# Split up chromosomes, with repaired and remaining breaks
	breakList.sort()
	chromFrags = splitChromosomes(chromosomes,breakList)
	chromList = copy.deepcopy(chromFrags)
	rings = []

	# Rejoin them
	for b1,b2 in repairs:
		c1, p1, d1, l1 = b1
		c2, p2, d2, l1 = b2
		n1, end1 = indexChrom(chromList,c1,p1,d1)
		n2, end2 = indexChrom(chromList,c2,p2,d2)

		# If ends are both the same chromosome, it's a ring
		if n1==n2:
			chromOne = chromList.pop(n1)
			rings.append(chromOne)
			continue

		# If not, stick them together and append new fragment
		if n1>n2:
			chromOne = chromList.pop(n1)
			chromTwo = chromList.pop(n2)
		else:
			chromTwo = chromList.pop(n2)
			chromOne = chromList.pop(n1)
		newFrag = appendFragments(chromOne,chromTwo,end1,end2)
		chromList.append(newFrag)

	linearChromosomes = characteriseChroms(chromList)
	ringChromosomes = characteriseChroms(rings,True)

	# Plot and print stats if requested
	if plot: plotAberrations.drawChroms(chromosomes, linearChromosomes, ringChromosomes, 
								   		inFile=inFile, outFile=outFile)

	print(index, '\t', breaks, '\t', len(remBreaks), 
		  '\t', misrepairStats(repairs, chromosomes),
	      '\t', calculateComplexities(chromList+rings), 
	      '\t', centricCheck(chromosomes,linearChromosomes),
	      '\t', centricCheck(chromosomes,ringChromosomes), 
	      end='')
	if baseBreaks!=None:
		lostDNA, lostFragments = dnaLoss(chromosomes, linearChromosomes+ringChromosomes)
		print('\t',fragmentDistribution(chromosomes,baseBreaks),
			  '\t',lostDNA, end='')
		if allFragments:
			print('\t\t','\t'.join(map(str,sorted(lostFragments))),end='')
		
	print()

	return chromList,rings, lostFragments


#######################################################
#
# Code to handle import from file in standard way
#
#######################################################
def repairFromFile(dataFile,doPlot=True):
	repairSets = []
	with open(dataFile) as dataSet:
		# First row is chromosome sizes. 
		# Chromosome defined as [ChromID, Start, End]
		# Measured in MBP, although this is normalised for plot
		chromData = next(dataSet)
		while chromData[0]=="#":
			chromData = next(dataSet)
		chromData = chromData.strip().split(",")
		chromData = [ [n,0,float(c)] for n,c in enumerate(chromData)]

		# Remaining rows are repair events. Defined as a series of repair pairs.
		# Blank lines denote separate exposures.
		currentSet = []
		for row in dataSet:
			if row[0]=="#":
				continue
			if len(row.strip())==0:
				if len(currentSet)>0:
					repairSets.append(currentSet)
				currentSet = []
				continue

			rowData = row.strip().rstrip(",").split(",")
			rowData = list(map(float,rowData))
			rowData = [ [rowData[0],rowData[1],rowData[2]],
			            [rowData[3],rowData[4],rowData[5]] ]
			currentSet.append(rowData)

		if len(currentSet)>0:
			repairSets.append(currentSet)

	# Finally, plot data
	for n,rep in enumerate(repairSets):
		fileIn = dataFile+" "+str(n)
		fileOut = dataFile.split('.')[0]+" "+str(n)+" aberrations.png"
		leftChroms, rings = doRepair(chromData,rep, plot=doPlot, inFile = fileIn, outFile=fileOut)
