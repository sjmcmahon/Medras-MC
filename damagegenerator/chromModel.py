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
import math
import random

chromCentres = []
radius = 1.0

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
	
# Subdivide a plane between minZ and maxZ into chroms chromosome points
def subDividePlane(minZ,maxZ,planeChroms):
	jitter = 0.05
	midZ = 0.5*(minZ+maxZ)
	zSpread = maxZ-minZ
	scaledR = np.sqrt(1-midZ*midZ)

	# Split into a number of rows, approximately equal to square root of chromosomes
	chromRows = int(pow(planeChroms,0.5))
	if random.random()>0.5:
		chromRows+=1

	# Assign at least one chromosome in each row
	rowSizes = [1]*chromRows
	for n in range(planeChroms-chromRows):
		rowSizes[random.randint(0,chromRows-1)]+=1
	
	# Calculate final positions for each chromosome
	cumChroms = 0
	currY = 0
	lastVol = 0
	positions = []
	for rowChroms in rowSizes:
		# Again, fraction of area is equal to by fraction of chromosomes
		# Split along Y dimension here.
		cumChroms += rowChroms
		newArea = cumChroms/planeChroms*(np.pi)

		# Guess initial Y value, then use root-finding to solve.
		# No simple polynomial, so use more complex method
		yGuess = 2*np.sqrt(newArea/np.pi)-1
		newY = (scipy.optimize.broyden1(lambda x:areaFunc(newArea,x),yGuess) ).item()

		# Calculate midpoint of Y slice, and X and Y spans
		midY = 0.5*(currY+newY)
		ySpan = newY-currY
		xSpan = 2*np.sqrt(1-midY*midY)

		# Distribute chromosomes along X dimension
		xPos = [( (x+1)/(rowChroms+1.0)-0.5)*xSpan for x in range(rowChroms)]

		# Generate new positions. Add small jitter to points to break regularity.
		positions=positions+ [ [scaledR*(x+jitter*xSpan/(rowChroms+1.0)*random.uniform(-0.5,0.5) ),
								scaledR*(midY+jitter*ySpan*random.uniform(-0.5,0.5)),
								midZ+zSpread*jitter*random.uniform(-0.5,0.5)] for x in xPos]

		currY = newY
	return positions

# Sub-divide a sphere into a number of chromosome territories
# Very simplified, no strong biological rationale
def subDivideSphere(noChrom,newRadius = 1.0):
	global chromCentres
	global radius
	radius = newRadius

	# Divide into a number of Z chromosome 'planes', scaled as cube root of noChrom. 
	# Floor or ceiling this value randomly.
	chromPlanes = int(pow(noChrom,1/3.0))
	if random.random()>0.5:
		chromPlanes+=1

	# Place at least 1 chromosome in each plane, then randomly distribute the rest
	chromPerPlane = [1]*chromPlanes
	for n in range(noChrom-chromPlanes):
		chromPerPlane[random.randint(0,chromPlanes-1)]+=1
	
	# Assign thickness of each Z plane based on chromosome number
	cumChroms = 0
	currZ = 0
	lastVol = 0
	planeZ = []
	for chroms in chromPerPlane:
		# Fraction of volume is given by current chromosome count over total chromosomes
		cumChroms += chroms
		newVol = cumChroms/noChrom*(4.0/3.0*np.pi)

		# Solve volume equation to get Z position corresponding to plane
		coeffs = [-np.pi/3, np.pi, 0, 2*np.pi/3-newVol]
		newZ = np.roots(coeffs)[1].real

		# Store starting and final Z, and chromosomes per plane
		planeZ.append([currZ,newZ,chroms])
		currZ = newZ

	# Distribute chromosomes in each plane
	chromCentres = []
	for row in planeZ:
		chromCentres=chromCentres+subDividePlane(*row)

	# Shuffle order and rotate to randomize chromosome positions
	random.shuffle(chromCentres)
	chromCentres = applyRotation(chromCentres)

	# Scale to final radius
	chromCentres = [np.array(c)*radius for c in chromCentres]

# Generate position of DNA, assuming DNA content linearly increases with increasing Z. 
# We don't know exact chromosome extent, so approximate this and fuzz on edges
def generateDNAPosition(x,y,z,c):
	chromCount = len(chromCentres)
	effRadius = radius/pow(chromCount,1.0/3.0)
	chromZ = c[2]
	deltaZ = z-chromZ
	if abs(deltaZ)<0.95*effRadius:
		return 0.5+(effRadius*effRadius*deltaZ-pow(deltaZ,3)/3.0)/(4.0/3.0*pow(effRadius,3.0))

	if deltaZ>0:
		return 1-np.exp(-6.6526/effRadius*deltaZ)

	if deltaZ<0:
		return np.exp(-6.6526/effRadius*abs(deltaZ))

# Generate chromosome data for a random x,y,z position. Simple nearest neighbour model.
# Assume this is in G1, so always chromatid 1. 
def modelChromosome(x,y,z):
	# Iterate through chromosomes to find nearest
	nearC = -1
	nearDist = 1E9
	for n, p in enumerate(chromCentres):
		distSq = ( (x-p[0])*(x-p[0])+ 
				   (y-p[1])*(y-p[1])+ 
				   (z-p[2])*(z-p[2]) )
		if distSq<nearDist:
			nearDist=distSq
			nearC = n

	# Return values here. Calculate DNA position based on Z position of event.
	chromosomeID = nearC
	chromatidID = 1
	chromosomeFrac = generateDNAPosition(x,y,z,chromCentres[nearC])
	chromosomeArm = 0 if chromosomeFrac<0.33 else 1 # Based on average centromere position

	return ["0, "+str(chromosomeID)+','+str(chromatidID)+','+str(chromosomeArm), chromosomeFrac]
