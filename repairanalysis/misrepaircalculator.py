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
import copy

import time

# Input parameters
# Physical DSB interaction rates, based on break complexity (per hour)
fastRate = 2.07
slowRate = 0.259

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
	seps = scipy.spatial.distance.pdist(positions)
	seps = scipy.spatial.distance.squareform(seps)
	rateTable = np.exp(-(seps*seps)/sigmaSq)
	rateTable = np.clip(rateTable, 1E-200, None)
	np.fill_diagonal(rateTable,0)

	return rateTable

# Pick a repair based on random sampling from weights
def pickRepair(interactionArray,n):
	cumuInteraction = np.cumsum(interactionArray)
	interactionSample = np.random.uniform()*cumuInteraction[-1]
	offsetInteraction = cumuInteraction-(interactionSample)
	chosenInteraction = np.where(offsetInteraction>=0)[0][0]

	return chosenInteraction

def singleRepair(breakList,rateTable,sigma=None, finalTime = np.inf):
	# Sort breaks by order of creation in time and set up interaction rates if needed
	if rateTable is None: 
		breakList.sort(key = lambda x:x[7])
		rateTable=buildRateTable(breakList,sigma)

	# Get fast/slow kinetic data from breaklist
	repairRate = np.array([fastRate/2 if b[2]==0 else slowRate/2 for b in breakList])

	# Sample interaction time for every break
	interactionSamples = -np.log(np.random.rand(len(breakList)))

	# Get base time, and initialise lists of live and pending breaks
	baseTime = breakList[0][7]
	liveBreaks = [n for n,b in enumerate(breakList) if b[7]<=baseTime]
	pendingBreaks = [n for n,b in enumerate(breakList) if b[7]>baseTime]
	lastBreak = liveBreaks[-1]
	if len(pendingBreaks)>0:
		nextBreakTime = breakList[pendingBreaks[0]][7]
	else:
		nextBreakTime = np.inf

	# Reporting variables
	repairEvents = []
	misrepairedPairs = []

	# Simulate repair until no more breaks remain, or we exceed the simulation limit
	while len(liveBreaks)+len(pendingBreaks)>0:
		if len(liveBreaks)>0:
			# For every break, baseline rate is repair rate divided by total interaction with all DSB ends.
			rateSums = np.sum(rateTable[:,0:lastBreak+1],axis=1)
			interactionTimes = np.array(interactionSamples)[liveBreaks]/(repairRate[liveBreaks]*rateSums[liveBreaks])
			nextTime = baseTime+min(interactionTimes)
		else:
			nextTime = nextBreakTime

		# Break loop if next event is after the simulation end
		if min(nextTime, nextBreakTime)> finalTime: break

		# If next repair is before next break, log the repair and remove the break ends
		if nextTime < nextBreakTime:
			# Get repaired break ends based on matching repair time, then select a pair
			endOne = liveBreaks[np.argmin(interactionTimes)]
			endTwo = pickRepair(rateTable[endOne,0:lastBreak+1],endOne)

			# Assign complexity for foci clearance based on most complex end
			complexity = max(breakList[endOne][2],breakList[endTwo][2])

			# Identify type of repair, and log details
			if breakList[endOne][0]!=breakList[endTwo][0]:
				# If ends aren't from the same break, then it's a misrepair. Log details. 
				if (breakList[endOne][3][1]==breakList[endTwo][3][1] and 
					breakList[endOne][3][2]==breakList[endTwo][3][2]):
					interChrom = 0
				else:
					interChrom = 1
				separation = np.sqrt(distanceToSq(breakList[endOne][1], breakList[endTwo][1]))
				misrepairedPairs.append([breakList[endOne], breakList[endTwo], separation, interChrom])

 				# Track if either end is the second end from a DSB. Log reduction in breaks if so
				p1Partner = endOne +1 - 2*(endOne%2)
				if breakList[p1Partner] == 0: repairEvents.append([nextTime,endOne,endTwo,complexity])
				p2Partner = endTwo +1 - 2*(endTwo%2)
				if breakList[p2Partner] == 0: repairEvents.append([nextTime,endOne,endTwo,complexity])

			else:
				# Matched DSB ends always clear 1 DSB, no misrepair
				repairEvents.append([nextTime,endOne,endTwo,complexity])

			# Tidy up break data
			liveBreaks.pop(liveBreaks.index(endOne))
			liveBreaks.pop(liveBreaks.index(endTwo))
			rateTable[:,endOne]=0
			rateTable[:,endTwo]=0
			breakList[endOne]=0
			breakList[endTwo]=0
		else:
			# If next DSB is before next repair, add pending DSBs and update interactions
			baseTime = nextBreakTime
			liveBreaks.append(pendingBreaks.pop(0))

			while len(pendingBreaks)>0 and baseTime>=breakList[pendingBreaks[0]][7]:
				liveBreaks.append(pendingBreaks.pop(0))

			interactionSamples = -np.log(np.random.rand(len(breakList)))
			if len(pendingBreaks)>0:
				nextBreakTime = breakList[pendingBreaks[0]][7]
			else:
				nextBreakTime = np.inf
			lastBreak = liveBreaks[-1]

	# Return a list of unrepaired breaks, if we halted before all were repaired
	if finalTime<np.inf:
		remBreaks = [b for b in breakList if b!=0]
		return misrepairedPairs, repairEvents, remBreaks

	# Otherwise, clean up any breaks which were missed and return 'full' repair
	# This is an effectively random pairing, but should hopefully just be tidying up corner cases
	if len(liveBreaks)+len(pendingBreaks)>0:
		remBreaks = [p for p in range(len(breakList)) if breakList[p]!=0]
		while len(remBreaks)>0:
			p1 = remBreaks.pop()
			p2 = remBreaks.pop()
			repairEvents.append([1E6,p1,p2,1])
			misrepairedPairs.append([breakList[p1],breakList[p2],-1,0])

	# Otherwise just return 'full' repair
	return misrepairedPairs, repairEvents, []

# Full repair in single pass.
def fullRepair(baseBreaks, sigma, repeats=1, addFociClearance=True, radius=1, 
			   chromSizes=None, sizeLimit=0, finalTime=np.inf):
	# Sort breaks in order of time of creation and build table of interaction rates
	baseBreaks.sort(key = lambda x:x[7])
	rateTable = buildRateTable(baseBreaks,sigma)

	fullMisrepairPairs = []
	fullRepairEvents = []

	for n in range(repeats):
		pairRates = []
		breakList = copy.deepcopy(baseBreaks)
		misrepairedPairs, repairEvents, remBreaks = singleRepair(breakList,rateTable.copy())
		
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
			fullRepairTimes.append(t-np.log(np.random.random())/repRate)

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
