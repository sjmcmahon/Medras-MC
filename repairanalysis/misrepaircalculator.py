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
# Core repair model code, which simulates repairs of initial 
# break distributions
#
# ############################################################################

import numpy as np
import scipy
import random
import copy

import time

# Input parameters
# Physical DSB interaction rates, based on break complexity (per hour)
fastRate = 2.07
slowRate = 0.259
relSlow = slowRate/fastRate

# Foci clearance rates (per hour)
fastFoci = 8.12
slowFoci = 0.405
mmejFoci = 0.0176

#Square of distance to
def distanceToSq(p1,p2):
	dx=p1[0]-p2[0]
	dy=p1[1]-p2[1]
	dz=p1[2]-p2[2]

	return dx*dx+dy*dy+dz*dz

# Interaction rate of two points
def interactionRate(p1,p2,sigma):
	return np.exp(logInteractionRate(p1,p2,sigma))

# Log interaction rate.
def logInteractionRate(p1,p2,sigma):
	return -distanceToSq(p1,p2)/(2*sigma*sigma)

def buildRateTable(baseBreaks,sigma):
	sigmaSq = 2*sigma*sigma
	positions = np.array([b[1] for b in baseBreaks])
	distances = scipy.spatial.distance.cdist(positions,positions)
	rateTable = np.maximum(np.exp(-np.power(distances,2)/sigmaSq),1E-200)
	positions = None
	distances = None
	complexity = np.array([[max(a[2],b[2]) for b in baseBreaks] for a in baseBreaks])
	rateTable = rateTable*np.power(relSlow,complexity)
	complexity = None
	return rateTable

def singleRepair(breakList,rateTable,sigma=None, finalTime = np.inf):
	if rateTable is None: rateTable=buildRateTable(breakList,sigma)

	# Calculate MC interaction rates, as -ln(rand())/interactionRate
	pairRates = [[i,j,-np.log(random.random())/rateTable[i][j] ] 
				 for i in range(len(breakList)) for j in range(i+1,len(breakList)) 
				 if rateTable[i][j]>1E-6 ]

	# Sort backwards, pop from end for efficiency
	pairRates.sort(key=lambda x: -x[2])
	totalBreaks = len(breakList)*0.5
	totalDSBs   = len(breakList)*0.5
	misrepBreaks = 0.0

	repairEvents = []
	misrepairedPairs = []

	# For each sorted pair, if they still exist, record the repair
	while totalBreaks>0 and len(pairRates)>0:
		p1,p2,t = pairRates.pop()
		if t>finalTime: break

		if breakList[p1]!=0 and breakList[p2]!=0:
			complexity = max(breakList[p1][2],breakList[p2][2])
			if breakList[p1][0]!=breakList[p2][0]:
				misrepBreaks+=2
				if (breakList[p1][3][1]==breakList[p2][3][1] and 
					breakList[p1][3][2]==breakList[p2][3][2]):
					interChrom = 0
				else:
					interChrom = 1
				separation = np.sqrt(distanceToSq(breakList[p1][1], breakList[p2][1]))
				misrepairedPairs.append([breakList[p1], breakList[p2], separation, interChrom])

 				# Track if either of these are the second repair from a DSB. Log if so
				p1Partner = p1 +1 - 2*(p1%2)
				if breakList[p1Partner] == 0:
					totalDSBs-=1
					repairEvents.append([t/fastRate,p1,p2,complexity])

				p2Partner = p2 +1 - 2*(p2%2)
				if breakList[p2Partner] == 0:
					totalDSBs-=1			
					repairEvents.append([t/fastRate,p1,p2,complexity])

			else:
				# Matched DSB ends always clear 1 DSB
				totalDSBs-=1
				repairEvents.append([t/fastRate,p1,p2,complexity])

			# Tidy up breaklist
			breakList[p1]=0
			breakList[p2]=0
			totalBreaks-=1

	# Return a list of unrepaired breaks, if we aborted before all were repaired
	if finalTime<np.inf:
		remBreaks = [b for b in breakList if b!=0]
		return misrepairedPairs, repairEvents, remBreaks

	# Otherwise, clean up any breaks which were missed and return 'full' repair
	if totalBreaks>0:
		remBreaks = [p for p in range(len(breakList)) if breakList[p]!=0]
		while len(remBreaks)>0:
			p1 = remBreaks.pop()
			p2 = remBreaks.pop()
			repairEvents.append([1E6,p1,p2,1])
			misrepairedPairs.append([breakList[p1],breakList[p2],-1,0])

	# Otherwise just return 'full' repair
	return misrepairedPairs, repairEvents

# Full repair in single pass.
def fullRepair(baseBreaks, sigma, repeats=1, addFociClearance=True, radius=1, 
			   chromSizes=None, sizeLimit=0, finalTime=np.inf):
	# Build table of interaction rates
	rateTable = buildRateTable(baseBreaks,sigma)

	fullMisrepairPairs = []
	fullRepairEvents = []

	for n in range(repeats):
		pairRates = []
		breakList = copy.deepcopy(baseBreaks)
		misrepairedPairs, repairEvents = singleRepair(breakList,rateTable)
		
		# Filter out any intra-chromosome breaks below size limit
		if sizeLimit>=0 and chromSizes is not None:
			for i in range(len(misrepairedPairs)-1,-1,-1):
				misrepair = misrepairedPairs[i]
				if misrepair[0][3][1]==misrepair[1][3][1]:
					chromID = misrepair[0][3][1]
					misrepSize = chromSizes[chromID]*abs(misrepair[0][4]-misrepair[1][4])
					if misrepSize<sizeLimit:
						misrepairedPairs.pop(i)

		fullMisrepairPairs.append(misrepairedPairs)
		fullRepairEvents+=repairEvents

	if addFociClearance:
		fullRepairTimes = []
		for t,p1,p2,complexity in fullRepairEvents:
			repRate = fastFoci
			if baseBreaks[p1][0]!=baseBreaks[p2][0]:
				repRate=mmejFoci
			else:
				if complexity>0:
					repRate = slowFoci
			fullRepairTimes.append(t-np.log(random.random())/repRate)

		fullRepairTimes = sorted(fullRepairTimes)
	else:
		fullRepairTimes = sorted([x[0] for x in fullRepairEvents])

	# Calculate rate of misrepair and stdev
	misrepairCounts = [len(m) for m in fullMisrepairPairs]
	misrepRate = np.mean(misrepairCounts)/(0.5*len(baseBreaks))
	stdevRate = np.std(misrepairCounts)/(0.5*len(baseBreaks))

	# Calculate inter-chromosome rate
	interChromEvents = sum(sum(misrep[3] for misrep in repeat) for repeat in fullMisrepairPairs)
	interChromRate = interChromEvents/max(1.0,1.0*sum(misrepairCounts))

	analyticRate, errRate = analyticRepair(baseBreaks,rateTable,sigma,radius)

	return misrepRate,stdevRate,interChromRate,fullRepairTimes, analyticRate, errRate

# Analytic repair method - gives fast approximation, but only valid for X-rays
def analyticRepair(breakList,rateTable,sigma,radius):
	correctRate = 0.0
	errRate = 0.0
	for i in range(len(breakList)):
		for j in range(i+1,len(breakList)):
			if breakList[i][0]!=breakList[j][0]:
				errRate += rateTable[i][j]

	errRate = errRate/(len(breakList)*0.5)

	base = 0.81026
	rate = 8.51
	skewCorrection = base+(1-base)*(1-np.exp(-rate*sigma/radius))
	errRate = errRate/skewCorrection

	if errRate<1E-10: return 0,errRate

	analyticMisrep = 1-2*(np.arctan(errRate+1)-np.arctan(1))/errRate
	return  analyticMisrep, errRate
