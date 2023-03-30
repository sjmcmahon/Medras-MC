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
# Initial version prepared by Stephen McMahon, implementing standard as 
# described in:
# 
# A New Standard DNA Damage (SDD) Data Format, Schuemann et al, Radiarion 
# Research, 2019
# https://doi.org/10.1667/RR15209.1
#
# Import and call sddparser.parseSDDFile(filename). The header, and each event
# are stored as dictionaries. Events are returned as a single object, which is
# a list of exposures, each of which is a list of individual damage events
#
# ############################################################################

import itertools
import re

# Read semicolon separated lines, stripping out anything between ## with regex
def delimitedRead(fileChain,delimiter = ';'):
    while True:
        line = ''.join(itertools.takewhile(lambda x: x!=delimiter,fileChain))
        if not line:
            break
        line = line.strip()
        line = re.sub(r'#.+?#', '', line)

        yield line.strip()
    
# Header and parsing helper functions for some particular formats
def parseEnergies(energyVals):
    energyList = (",".join(energyVals)).split('/')
    energyList = [e.split(',') for e in energyList]
    return energyList

def parseVolumes(volumeVals):
    worldVolume = [volumeVals[0]]+list(map(float,volumeVals[1:7]))
    scoringVolume = [volumeVals[7]]+list(map(float,volumeVals[8:]))
    return worldVolume,scoringVolume

def parseChromosomes(chromosomeVals):
    return [int(chromosomeVals[0]),list(map(float,chromosomeVals[1:]))]

def parsePosition(posString):
    posGroups = posString.split('/') # Split into up to three groups of coordinates 
                                     # (Centre and 2x bounding XYZ sets)
    posVals = [ list(map(float,s.split(','))) for s in posGroups] # Split each into XYZ
    posVals = posVals + [ [] ]*(3-len(posVals))                   # Pad as needed
    return posVals

def parseProliferation(prolifString):
    if len(prolifString)>1:
        return [int(prolifString[0]), prolifString[1] ]
    return [int(prolifString[0]), [] ]

# General array parser - parse value A using specificed function.
# We use try/except block to catach a few different possible line formats
# Throw "None" if something goes wrong somewhere
def parseGeneral(A,function):
    if len(A)==1 and A[0].strip()=='': return function(-1) # Catch empty fields, return -1.
    
    # Try to do direct casting for single values. Catch stray units or trailing characters.
    if len(A)==1: 
        try:
            return function(A[0])
        except:
            try:
                return function(A[0].strip().split()[0])
            except:
                print('Error in parsing:',A)
                return None

    # As above, but for lists.
    try: 
        return list(map(function,A))
    except:
        retList = []
        for element in A:
            try:
                retList.append(function(element))
            except:
                if len(element)!=0:
                	print('Error in parsing:',element,' from: ',A)
                retList.append(None)
        return retList

# Parse header, using iter_tools chain file object
def parseHeader(fc):
    header = {}
    rawHeader = []
    for line in delimitedRead(fc):
        line = line.rstrip()
        if len(line)==0 or line[0]!='#':
            rawHeader.append(line)
            if len(rawHeader)>27 or line=="***EndOfHeader***":
                break

    # Split header into csv-fields
    for n,row in enumerate(rawHeader):
        rawHeader[n]=list(map(str.strip,row.split(',')))

    # Simple string fields are included directly
    stringFields = [[0,"SDD Version"],[1,"Software"],[2,"Author"],[3,"Simulation details"],
    			    [4,"Source"],[11,"Irradiation Target"],[25,"Additional Information"]]

    for field,fieldName in stringFields:
        header[fieldName] = "".join(rawHeader[field][1:])

    # A number of fields are arrays of either ints or floats and can be parsed together
    numberFields = [[6,'Incident Particles',int],[7,'Mean Particle Energy',float],
                    [9,'Particle Fraction',float],[10,'Dose or Fluence',float],[15,'DNA Density',float],
                    [16,'Cell Cycle Phase',float],[17,'DNA Structure',int],[18,'In vitro/in vivo',int],
                    [20,'Microenvironment',float],[21,'Damage definition',float],
                    [23,'Damage and Primary Count',int],[24,'Data entries',int]]
    for field,fieldName,function in numberFields:
        header[fieldName] = parseGeneral(rawHeader[field][1:],function)

    # A small number of fields are single elements or special formats
    header['Source Type'] = int(rawHeader[5][1])
    header['Energy Distribution'] = parseEnergies(rawHeader[8][1:])
    header['World Volume'],header['Scoring Volume'] = parseVolumes(rawHeader[13][1:])
    header['Proliferation status'] = parseProliferation(rawHeader[19][1:])
    header['Chromosomes'] = parseChromosomes(rawHeader[14][1:])
    header['Time'] = float(rawHeader[22][1])

    return header

# Parse single event, extracting fields identified in header
def parseEvent(e,fields):
	# New event field is mandatory
    event = {}
    event['NewEvent']=int(e[0][0])
    if len(e[0])>1:
        event['EID'] =int(e[0].split(',')[1])

    # Now, parse remaining fields. 
    # As with header, several are parsed in similar ways   
    standardFields = [[2,'Chromosome ID',int], [3,'Chromosome Position',float],[4,'Cause',int],
    				  [5,'Damage Types',int],[7,'DNA Sequence',str]]
    for field,fieldName,function in standardFields:
        if fields[field]:
            targetField = sum(fields[0:field])
            event[fieldName] = parseGeneral(e[targetField][0:].split(','),function)

    # Now, some with special requirements
    if fields[1]:
        targetField = sum(fields[0:1])
        event['Pos']=parsePosition(e[targetField])
    if fields[6]:
        targetField = sum(fields[0:6])
        event['Damage Spec'] = [ [d.strip() for d in damage.split(',')] 
        						  for damage in e[targetField].split('/') if len(damage)>0]
    if fields[8]:
        targetField = sum(fields[0:8])
        event['Lesion Time'] = list(map(float,e[targetField].split('/')))

    # Final set of fields describe all of the particles involved in damage
    if fields[9]:
        targetField = sum(fields[0:9])
        event['Particle Types'] = list(map(int,e[targetField].split(',')))
    if fields[10]:
        targetField = sum(fields[0:10])
        event['Energies'] = list(map(float,e[targetField].split(',')))
    if fields[11]:
        targetField = sum(fields[0:11])
        event['Translation'] = [ [float(t) for t in translation.split('/')] 
        						  for translation in e[targetField].split(',')]
    if fields[12]:
        targetField = sum(fields[0:12])
        event['Direction'] = [ [float(t) for t in translation.split('/')] 
        						for translation in e[targetField].split(',')]
    if fields[13]:
        targetField = sum(fields[0:13])
        event['Particle Time'] = list(map(float,e[targetField].split(',')))

    return event

# Parse all events in data block
def parseDataBlock(fc,fields):
    events = []
    currEvent = []
    blockSize = sum(fields)
    for field in delimitedRead(fc):
        currEvent.append(field)
        if len(currEvent) == blockSize:
            newData = parseEvent(currEvent,fields)
            if newData['NewEvent']==2:
                events.append([])
            events[-1].append(newData)
            currEvent = []

    return events

# General parsing method, retuns header as dictionary, and list of events as dictionaries
def parseSDDFile(fileName,verbose=False):
    with open(fileName) as f:
        fc = itertools.chain.from_iterable(f)
        header = parseHeader(fc)
        fields = header['Data entries']
        events = parseDataBlock(fc,fields)

        if (verbose and header['Damage and Primary Count']!=-1 and 
            sum([len(e) for e in events])!=header['Damage and Primary Count'][0]):
            print('Event count mismatch. Expected: ', header['Damage and Primary Count'], '; read: ',len(events))
        return header,events
