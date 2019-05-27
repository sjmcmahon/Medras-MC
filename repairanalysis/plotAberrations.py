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

import os
import matplotlib as mpl
haveDisplay = "DISPLAY" in os.environ
if not haveDisplay:
	mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches
import numpy as np

#######################################################
#
# Methods to prepare MFish Figure
#
#######################################################
def chromLength(chromosome):
	return sum(abs(f[2]-f[1]) for f in chromosome)

# Count centromeres in a single fragment
# Centromere at a fractional position of pos
def centricCount(c, baseChromosomes, pos=0.5):
	centromereCount=0
	for f in c[3]:
		start = min(f[1],f[2])
		end = max(f[1],f[2])
		cent = pos*baseChromosomes[f[0]][2]
		if cent>=start and cent<=end:
			centromereCount+=1
	return centromereCount

# Count centromeres in each chromosome and return appropriate tag
def centromereTag(chromosome,baseChromosomes,pos=0.5):
	centromereCount = centricCount(chromosome,baseChromosomes,pos)

	if centromereCount==1:
		return ""
	if centromereCount==0:
		return "#"
	return "*"

# Colourset extracted from a randomly chosen MFish image, for 23 human chromosomes
def fetchColor(c):
	colorSet = [ (0.847058824,0.882352941,0.305882353),(0.749019608,0.007843137,0.211764706),(0.517647059,0.682352941,0.839215686),
	 			  (0.019607843,0.670588235,0.28627451),(0.68627451,0.262745098,0.57254902),(0.88627451,0.474509804,0.145098039),
	 			  (0.901960784,0.639215686,0.184313725),(0.48627451,0.431372549,0.596078431),(0.552941176,0.847058824,0.850980392),
	 			  (0.607843137,0.701960784,0.278431373),(0.168627451,0.654901961,0.545098039),(0.854901961,0.17254902,0.223529412),
	 			  (0.380392157,0.807843137,0.858823529),(0.949019608,0.08627451,0.184313725),(0.588235294,0.447058824,0.666666667),
	 			  (0.996078431,0.988235294,0.498039216),(0.270588235,0.650980392,0.831372549),(0.964705882,0.937254902,0.945098039),
	 			  (0.168627451,0.698039216,0.298039216),(0.301960784,0.435294118,0.674509804),(0.894117647,0.921568627,0.498039216),
	 			  (0.968627451,0.968627451,0.576470588),(0.462745098,0.788235294,0.301960784) ]
	return colorSet[c//2]

# Core drawing method
def drawChroms(baseChromosomes,chromosomes,rings, inFile=None, outFile=None):
	noRows = int(np.ceil(len(baseChromosomes)/10.0)) # 5 pairs per row
	chromSet = chromosomes+rings
	chromSet.sort(key=lambda x:x[1])

	plt.style.use('dark_background')
	ax = plt.axes()
	ax.set_facecolor('black')
	width = 10
	intraOffset = 5
	interOffset = 30
	currX = interOffset
	maxHeight = 0
	
	maxHeight = max(c[0] for c in chromSet)
	yOffset = maxHeight*1.1+10
	refY = -yOffset/2.0
	lastChrom = -1
	rowStart = -100
	if inFile is not None:
		ax.text(currX,refY-0.25*yOffset,inFile,color=(0.5,0.5,0.5))

	for c in chromSet:
		if lastChrom//2==c[1]//2:
			currX += intraOffset
		else:
			lastChrom=c[1]
			currX+=interOffset
			if c[1]-rowStart>10:
				refY -= yOffset
				currX = interOffset
				rowStart = c[1]

		currY=refY
		# Draw a linear chromosome
		if c[2]==False:
			totHeight = c[0]/2.0
			ax.text(currX+width/2.0, currY-yOffset*0.25,
					str((c[1]//2)+1)+centromereTag(c,baseChromosomes), ha="center",
					color=(0.5,0.5,0.5))
			for f in c[3]:
				fragHeight = abs(f[2]-f[1])/2.0
				rect = matplotlib.patches.Rectangle((currX,currY), width, fragHeight,
													color=fetchColor(f[0]))
				ax.add_patch(rect)
				currY += fragHeight
			currX += width
		else:
			# Draw a ring chromosome
			totLength = c[0]/2.0
			circRad = np.sqrt(totLength*width/3.14159)
			currX += max(circRad,width/2.0)
			ax.text(currX, currY-yOffset*.25,
				    str((c[1]//2)+1)+centromereTag(c,baseChromosomes), ha="center",
				    color=(0.5,0.5,0.5), fontweight="bold")
			currTheta = 0
			for f in c[3]:
				fracLength = abs(f[2]-f[1])/2.0
				deltaTheta = 360*(fracLength/totLength)
				wedge = matplotlib.patches.Wedge((currX,currY+circRad), circRad, currTheta,
												  currTheta+deltaTheta, color=fetchColor(f[0]))
				ax.add_patch(wedge)
				currTheta=currTheta+deltaTheta
			currX+=max(circRad,width/2.0)

	ax.autoscale()
	ax.axis('off')
	if outFile is None:
		plt.savefig("ChromAberrs.png")
	else:
		plt.savefig(outFile)

	if haveDisplay:
		plt.show()
	else:
		plt.close()