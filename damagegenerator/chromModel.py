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
import numpy as np
import scipy.optimize
import scipy.spatial
import math
import random

chromCentres = []
radius = 1.0
bpUncertainty = 0 # Variation in base pair position around a grid point

# Calculate area fraction associated with some x position for root-finding.
# Return illogical values if outside range to ensure good fit. 
def areaFunc(A,x):
	if x<-1: return -A-x/10
	if x>1:  return A+(x-1)/10

	return (x*np.sqrt(1-x*x)+np.arcsin(x)+np.pi/2) - A

# Return the rotation matrix associated with counterclockwise rotation about
# the given axis by theta radians.
def rotation_matrix(axis, theta):
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

# Generate a random rotation matrix to represent an arbitrary rotation around theta and phi
def applyRotation(points):
	# Generate random rotation vector
	theta = np.arccos(2*random.random()-1)
	phi = 2*np.pi*random.random()
	unitVector = [np.cos(phi)*np.sin(theta),np.sin(phi)*np.sin(theta),np.cos(theta)]

	# Apply random rotation matrix around this vector
	rot = 2*np.pi*random.random()
	rotMatrix = rotation_matrix(unitVector,rot)

	# Apply matrix to points and return
	return [np.dot(rotMatrix,p) for p in points]
	
# Sub-divide a sphere into a number of chromosome territories
# Very simplified, no strong biological rationale

# Assign chromosomes in a rigid grid, 4 chromosomes deep on each side of nucleus
# Can probably be generalised, but some tricky aspects
def assignChromosome(point,radius):
	# Sort by Z, thresholds at 0.36, 0.61, 0.8 and 1 radii. Calculate h based on this.
	x,y,z = [c/radius for c in point]
	if z<0:
		h = 1+z # Radius Minus (minus z)
		zOffset = 0
	else:
		h = 1-z
		zOffset = 23 # Z Offset for far side of nucleus

	if h<=0.36: # 4 chromosomes
		if y<0: 
			yOffset = 0 
		else:
			yOffset = 2
		if x<0: return 0 + zOffset + yOffset
		if x>=0: return 1 + zOffset + yOffset

	if h<=0.61: # 6 chromosomes
		if y<0: 
			yOffset = 0  
		else:
			yOffset = 3
		edge = 0.22 # Need to mark position with 1/3 on each side
		if x<=-edge: return 4 + zOffset + yOffset
		if x>-edge and x<edge: return 5 + zOffset + yOffset
		if x>=edge: return 6 + zOffset + yOffset

	if h<=0.8: # 6 chromosomes, but split along X first instead this time
		if x<0: 
			xOffset = 0  
		else:
			xOffset = 3
		edge = 0.245 # Need to mark position with 1/3 on each side
		if y<=-edge: return 10 + zOffset + xOffset
		if y>-edge and y<edge: return 11 + zOffset + xOffset
		if y>=edge: return 12 + zOffset + xOffset

	if h<=1: # 7 chromosomes, 4 below bound 3 above
		if y<0.14: 
			edge = 0.42
			if x<=-edge: return 16 + zOffset
			if x>-edge and x<=0: return 17 + zOffset 
			if x>0 and x<=edge: return 18 + zOffset
			if x>=edge: return 19 + zOffset
		else:
			edge = 0.25 # Need to mark position with 1/3 on each side
			if x<=-edge: return 20 + zOffset
			if x>-edge and x<edge: return 21 + zOffset
			if x>=edge: return 22 + zOffset 

	print('Failed to assign chromosome!',point)
	return -1

# Very crude 3D splitting to give regions of chromosome. Partition along X, Y, Z sequentially.
def partitionChromosome(chromPoints, depth=1):
	finalSplits = []

	xSort = chromPoints[np.argsort(chromPoints[:,0],0)]
	xPartition = [xSort[:len(xSort)//2],xSort[len(xSort)//2:]]
	for xSplit in xPartition:
		ySort = xSplit[np.argsort(xSplit[:,1],0)]
		yPartition = [ySort[:len(ySort)//2],ySort[len(ySort)//2:]]
		for ySplit in yPartition:
			zSort = ySplit[np.argsort(ySplit[:,2],0)]
			zPartition = [zSort[:len(zSort)//2],zSort[len(zSort)//2:]]

			for zSplit in zPartition:
				zSplit = zSplit[np.lexsort((zSplit[:,0],zSplit[:,1],zSplit[:,2]))]
				if depth>1:
					finalSplits.append(partitionChromosome(zSplit,depth-1))
				else:
					finalSplits.append(zSplit)
	
	return np.concatenate(finalSplits)

# Note - we don't use noChrom for now, held as placeholder for later
def subDivideSphere(noChrom = 46,newradius = 4.65, gridSize = 106):
	global chromCentres
	global chromTree
	global radius
	global bpUncertainty

	if noChrom != 46: print("Warning: Currently restricted to human (46 chromosome) nuclei")

	radius = newradius
	gridPoints = [(n-0.5*gridSize+0.5)*2*radius/gridSize for n in range(gridSize)]

	allC = []
	pointGrid = []
	for z in gridPoints:
		for x in gridPoints:
			for y in gridPoints:
				if np.linalg.norm([x,y,z])<=radius:
					c=assignChromosome([x,y,z],radius)
					allC.append(c)
					pointGrid.append([x,y,z,c])

	newPointGrid = []
	# Update with chromosome positions, just incrementing in order. Guarantees some locality.
	for c in range(46):
		chromPoints = np.array([np.array(p) for p in pointGrid if p[3]==c])
		chromPoints = partitionChromosome(chromPoints, depth=2)
		chromPosition = [(c+0.5)/len(chromPoints) for c in range(len(chromPoints))]
		chromPoints= np.c_[chromPoints, chromPosition]
		newPointGrid.append(chromPoints)

	# Sort and tidy everything
	newPointGrid = np.concatenate(newPointGrid)
	chromCentres = newPointGrid[np.lexsort((newPointGrid[:,0],newPointGrid[:,1],newPointGrid[:,2]))]

	# Apply a rotation
	rotPoints = applyRotation(chromCentres[:,0:3])
	chromCentres = np.c_[rotPoints, chromCentres[:,3:]]

	# Calculate average points per chromosome to estimate uncertainty
	bpUncertainty = 1/(len(chromCentres)/46)

	# Finally, build a KD tree to enable quick lookup of points
	chromTree = scipy.spatial.KDTree(chromCentres[:,0:3])

# Generate chromosome data for a random x,y,z position. Simple nearest neighbour model.
# Assume this is in G1, so always chromatid 1. 
def modelChromosome(x,y,z):
	# Iterate through chromosomes to find nearest
	nearDist = 1E9

	# Use the KD tree to quickly look up nearest chromosome
	d, idx = chromTree.query([x,y,z])
	chromosomeID = int(chromCentres[idx][3])
	chromosomeFrac = chromCentres[idx][4] + 0.5 * bpUncertainty * (1-np.random.uniform())

	# Return values here. Calculate DNA position based on Z position of event.
	chromatidID = 1
	chromosomeArm = 0 if chromosomeFrac<0.33 else 1 # Based on average centromere position

	return ["0, "+str(chromosomeID)+','+str(chromatidID)+','+str(chromosomeArm), chromosomeFrac]
