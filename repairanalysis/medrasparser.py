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
# Spectra parser for SDD v1.0. 
# Converts SDD-parsed data into MEDRAS format
# 
# A New Standard DNA Damage (SDD) Data Format, Schuemann et al, Radiarion 
# Research, 2019
# https://doi.org/10.1667/RR15209.1
#
# ############################################################################

import numpy as np
import random
import scipy.stats

from . import sddparser

# Make single point, in sphere 
def makePoint(r):
    x=r*(1-2*random.random())
    y=r*(1-2*random.random())
    z=r*(1-2*random.random())

    while x*x+y*y+z*z>r*r:
        x=r*(1-2*random.random())
        y=r*(1-2*random.random())
        z=r*(1-2*random.random())

    return np.array([x,y,z])

def separateDSBs(dsbList):
    prevPos = []
    for n in range(len(dsbList)-1,0,-2):
        testBreak = dsbList[n]
        pos = [testBreak[3][1], testBreak[4]]
        if pos in prevPos:
            offset = 1E-9*random.random()*testBreak[4]
            dsbList[n][4]+=offset
            dsbList[n-1][4]+=offset
        prevPos.append(pos)

# Parse spec file to break distribution for repair model
def parseToBreaks(fileName,sig,r_sep = 0.0001,basicStats=False,verbose=False):
    header,events = sddparser.parseSDDFile(fileName,verbose)
    breaks = [ ]
    complexities = []
    totalBreaks = 0
    emptySets = 0

    totalRows = 0
    totalBD=0
    totalSSB=0
    totalDSB=0
    simpleDSB=0
    complexDSB=0
    multiDSB=0

    maxChrom = 0

    for dataSet in events:
        complexity = 0
        breaks.append([])
        index = 0
        setDSB = 0
        newEventStatus = 0
        for e in dataSet:
            totalRows+=1
            totalBD+=e['Damage Types'][0]
            totalSSB+=e['Damage Types'][1]
            totalDSB+=e['Damage Types'][2]

            if e['NewEvent']>newEventStatus:
                newEventStatus=e['NewEvent']
            if e['Damage Types'][2]==1:
                setDSB+=1
                totalBreaks+=1

                pos = np.array(e['Pos'][0])
                if len(pos)>3: pos = pos[0:3]
                
                pos2 = pos + makePoint(r_sep)
                if sum(e['Damage Types'][0:2])>2:
                    complexBreak = 1
                    complexity +=1
                else:
                    complexBreak = 0
                
                if e['Damage Types'][2]>1:
                    multiDSB+=1
                else:
                    if complexBreak==1:
                        complexDSB+=1
                    else:
                        simpleDSB+=1

                chromID = [-1,-1,-1,-1]
                chromPos = 0
                if 'Chromosome ID' in e:
                    chromID = e['Chromosome ID']
                    if chromID[1]>maxChrom:
                        maxChrom=chromID[1]
                if 'Chromosome Position' in e:
                    chromPos = e['Chromosome Position']

                # Calculate damage time, cast to hours. Set individual times to 0 if not provided.
                lesionTime = min(e['Lesion Time']) if 'Lesion Time' in e else 0
                particleTime = min(e['Particle Time']) if 'Particle Time' in e else 0
                damageTime = (lesionTime+particleTime)/(60*60*1E9)

                if 'Cause' in e: 
                    cause = e['Cause']
                else:
                    cause = 0

                # Break is: Index, position, complexity, chromosome ID, upstream/downstream, new event status, time and cause.
                breaks[-1]+=[[index,pos,complexBreak,chromID[:],chromPos,-1,newEventStatus, damageTime, cause],
                             [index,pos2,complexBreak,chromID[:],chromPos,1,0, damageTime, cause]]
                index+=1
                newEventStatus=0
        complexities.append(complexity)
        # If we have zero or one DSBs, throw this away as there can be no binary misrepair
        if setDSB<1:
            breaks.pop()
            complexities.pop()
            emptySets+=1
        else:
            # Separate out DSBs whose genetic position overlaps
           separateDSBs(breaks[-1])

    # Get effective sigma, and dose if provided
    scaledSigma = sig*scipy.stats.gmean(header['Scoring Volume'][1:4])
    if header['Dose or Fluence'][0]==1:
        dose = header['Dose or Fluence'][1]
    else:
        dose = -1

    # If max chromosome is equal to number of chromosomes, we need to offset by 1 to correct indexing
    if maxChrom==int(header['Chromosomes'][0]):
        for breakList in breaks:
            for b in breakList:
                b[3][1]=b[3][1]-1

    complexFrac = [1.0*c/(len(b)/2.0) for c,b in zip(complexities,breaks)]
    totalComplex = sum(complexities)/(sum([len(b) for b in breaks])/2.0)

    if verbose:
        print(fileName.split('/')[-1],'\t', 'Read:', len(breaks), 
              ' sets of breaks, with average complexity ', np.mean(complexFrac), 
              ', stdev: ',  np.std(complexFrac), '. Also read ',emptySets, ' empty break sets.')
    return breaks,scaledSigma, header['Mean Particle Energy'], np.mean(complexFrac), np.std(complexFrac), dose, emptySets, header
