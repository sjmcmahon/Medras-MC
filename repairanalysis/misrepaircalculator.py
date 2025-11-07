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

# Repair failure parameters - derived from analytic Medras rates
repairFailRate = 1-0.98537
interChromRate = 0.1348
largeDelRate = 0.2426
largeDelRateIfIntra = largeDelRate/(1-interChromRate)

###################
# Utility functions
###################
# Square of distance to
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

# Get the partner end ID for a given other end ID
partnerBreak = lambda p: p+1 - 2*(p%2)

###################
# Analyze a set of misrepair events to generate kinetics
def calculateRepairKinetics(repairEvents, addFociClearance):
	repairTimes = []
	unrepairedEnds = []
	for t, p1, p2, complexity in repairEvents:
		breaksRepaired = 0

		# If this is a correct repair, clears one DSB, with complexity from event
		if partnerBreak(p1) == p2:
			breaksRepaired = 1
		else:
			# If it's a mismatch, can clear 0, 1 or 2 ends depending on other breaks
			# And set complexity to 2, to reflect MMEJ for mismatch
			complexity = 2
			for p in (p1, p2):
				p=abs(p)
				try:
					partner = partnerBreak(p)
					partnerIndex =  unrepairedEnds.index(partner) # Will except here if no match.
					breaksRepaired+=1
					unrepairedEnds.pop(partnerIndex)
				except:
					unrepairedEnds.append(p) # If no match, add to unrepaired list
					continue

		# For this number of breaks, append this time as repair time
		# optionally plus time for foci clearance
		for n in range(breaksRepaired):
			candidateTime = t
			if addFociClearance:
				repRate = fastFoci
				if min(p1,p2)<0 or complexity>1:
					repRate=mmejFoci
				else:
					if complexity>0:
						repRate = slowFoci
				candidateTime = t-np.log(np.random.random())/repRate
			repairTimes.append(candidateTime)

	return repairTimes

# Select new chromosome location for mistmatch event
def repairFailureLocation(currChrom, breakPos, chromSizes):
	largeDelSize = 3 # In MBP
	# Version 4 - Correct inter-chrom rate and intra-chrom rate in own function
	interChrom = 0
	if np.random.rand()<interChromRate:
		interChrom = 1
		# If inter-chromosome, put on a random point on randomly selected other chromosome
		chromID = np.random.choice([i for i in range(0, len(chromSizes)) if i != currChrom[1]])
		#chromID = np.random.randint(0,len(chromSizes))
		#while chromID == currChrom: chromID = np.random.randint(0,len(chromSizes)) # Protect against getting same chromosome
		newChrom = [0, chromID, np.random.randint(0,2), np.random.randint(0,2)]
		newBreakLoc = np.random.rand()
	else:
		# if intra-chromosome, test if it's a large or a small deletion and size appropriately
		newChrom = copy.copy(currChrom)
		chromSize = chromSizes[currChrom[1]]
		if chromSize<largeDelSize*2 or np.random.rand()>largeDelRateIfIntra:
			# If it's a small break
			newBreakLoc = np.random.uniform(max(0,breakPos-largeDelSize/chromSize),min(1,breakPos+largeDelSize/chromSize)) # Sample within 3 MBP
		else:
			# Sample randomly outside a 3 MBP window around the break
			newBreakLoc = np.random.uniform(0,(chromSize-2*largeDelSize)/chromSize) # Sample length of chromosome minus 2*Large del size window
			if newBreakLoc > breakPos-largeDelSize/chromSize: newBreakLoc += 2*largeDelSize/chromSize # Then offset if after break

	return newChrom, newBreakLoc, interChrom

# Add single-event misrepairs, which fail separately from binary events
def addRepairFailure(misrepairEvents, repairEvents, breakList, chromSizes):
	newBreaks = []
	newRepairs = []
	newMisreps = []

	nextBreak = len(breakList)+1
	popList = []
	for n, repEvent in enumerate(repairEvents):
		if partnerBreak(repEvent[1]) == repEvent[2]: # If this is a correct repair, sample for repair failure
			if np.random.rand()<repairFailRate:
				popList.append(n) 			# Track list of which breaks to remove at end
				endOne = breakList[repEvent[1]]
				endTwo = breakList[repEvent[2]]

				newChrom, newBreakLoc, interChrom = repairFailureLocation(endOne[3], endOne[4], chromSizes) # Sample location of new misrepair

				upDown = [1,-1] if np.random.rand()>0.5 else [-1,1] # Randomly set orientation of breaks relative to centromere

				# Add breaks - index, position, complexity, chromosome ID, chromosome pos, centromere side, new event status, time and cause
				newBreaks += [[-(nextBreak+1), endOne[1],endOne[2], newChrom, newBreakLoc, upDown[0], endOne[6], endOne[7], endOne[8]],
							  [-(nextBreak+2), endTwo[1],endTwo[2], newChrom, newBreakLoc, upDown[1], endTwo[6], endTwo[7], endTwo[8]] ]
				nextBreak += 2

				# Record repair event - time, end 1 ID, end 2 ID, complexity
				newRepairs += [[repEvent[0], repEvent[1], newBreaks[-2][0], endOne[2]],
							   [repEvent[0], repEvent[2], newBreaks[-1][0], endTwo[2]]]

				# Record misrepair event - two break ends, separation which we log as 0, and inter-chromosome flag
				newMisreps += [[endOne, newBreaks[-2], 0, interChrom],
							   [endTwo, newBreaks[-1], 0, interChrom] ]

	for n in sorted(popList, reverse=True): repairEvents.pop(n) # Pop repair events which are to be replaced

	return misrepairEvents+newMisreps, repairEvents+newRepairs, breakList+newBreaks

def singleRepair(initBreakList, rateTable, sigma=None, finalTime = np.inf, chromSizes = None, repairFailure = True):
	# Sort breaks by order of creation in time and set up interaction rates if needed
	initBreakList.sort(key = lambda x:x[7])
	breakList = copy.deepcopy(initBreakList)
	if rateTable is None:
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

			# If it's a misrepair, log it and appropriate details
			if breakList[endOne][0]!=breakList[endTwo][0]:
				# If ends aren't from the same break, then it's a misrepair. Log details.
				if (breakList[endOne][3][1]==breakList[endTwo][3][1] and
					breakList[endOne][3][2]==breakList[endTwo][3][2]):
					interChrom = 0
				else:
					interChrom = 1
				separation = np.sqrt(distanceToSq(breakList[endOne][1], breakList[endTwo][1]))
				misrepairedPairs.append([breakList[endOne], breakList[endTwo], separation, interChrom])

			# Record all repair events for future kinetics analysis
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

	# Once all binary misrepair is done, apply repair failure events if requested
	if repairFailure:
		misrepairedPairs, repairEvents, newBreakList = addRepairFailure(misrepairedPairs, repairEvents, initBreakList, chromSizes)

	# Return a list of unrepaired breaks, if we halted before all were repaired
	if finalTime<np.inf:
		remBreaks = [b for b in breakList if b!=0]
		return misrepairedPairs, repairEvents, remBreaks

	# If full repair is requested clean up any breaks which were missed
	# This is a random pairing, but should hopefully just be tidying up corner cases
	if len(liveBreaks)+len(pendingBreaks)>0:
		remBreaks = [p for p in range(len(breakList)) if breakList[p]!=0]
		while len(remBreaks)>0:
			p1 = remBreaks.pop()
			p2 = remBreaks.pop()
			repairEvents.append([1E6,p1,p2,1])
			misrepairedPairs.append([breakList[p1],breakList[p2],-1,0])

	# Return final data
	return misrepairedPairs, repairEvents, []

# Full repair in single pass.
def fullRepair(baseBreaks, sigma, repeats=1, addFociClearance=True, radius=1,
			   chromSizes=None, sizeLimit=0, finalTime=np.inf, repairFailure = True):
	# Sort breaks in order of time of creation and build table of interaction rates
	baseBreaks.sort(key = lambda x:x[7])
	rateTable = buildRateTable(baseBreaks,sigma)

	# Data stores
	fullMisrepairPairs = []
	fullRepairTimes = []

	# Iterate over repeats
	for n in range(repeats):
		pairRates = []
		breakList = copy.deepcopy(baseBreaks)
		misrepairedPairs, repairEvents, remBreaks = singleRepair(breakList,rateTable.copy(),
																 chromSizes = chromSizes, repairFailure = repairFailure)

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
		fullRepairTimes += calculateRepairKinetics(repairEvents,addFociClearance)

	# Calculate rate of misrepair and stdev
	misrepairCounts = [len(m) for m in fullMisrepairPairs]
	misrepRate = np.mean(misrepairCounts)/(0.5*len(baseBreaks))
	stdevRate = np.std(misrepairCounts)/(0.5*len(baseBreaks))

	# Calculate inter-chromosome rate
	interChromEvents = sum(sum(misrep[3] for misrep in repeat) for repeat in fullMisrepairPairs)
	interChromRate = interChromEvents/max(1.0,1.0*sum(misrepairCounts))

	analyticRate, errRate = analyticRepair(baseBreaks,rateTable,sigma,radius)

	return misrepRate,stdevRate,interChromRate, fullRepairTimes, analyticRate, errRate

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
