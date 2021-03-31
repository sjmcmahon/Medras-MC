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

# Write SDD header
def writeHeader(f, energy, writeSparse, incident, dose, geometry, DNADensity, bdRange, damageCount, primaryCount):
	sourceType = 1
	radius = geometry[1]
	f.write("SDD version, SDDv1.0;\n")
	f.write("Software, McMahon Empirical Model. v0.4;\n")
	f.write("Author, Stephen J McMahon, stephen.mcmahon@qub.ac.uk, 12/06/2018, McMahon et al, Scientific Reports 10790;\n")
	f.write("Simulation Details, Empirical model of DNA damage distribution based on radial energy. Nucleus radius "+str(radius)+" um giving DSB energy of 60.1 keV/DSB on average;\n")
	f.write("Source, Plane-parallel particle beam uniformly exposing nucleus. This may be composed of one more particle types, identified below;\n")
	f.write("Source type, "+str(sourceType)+";\n")
	if type(incident) is list: incident = ','.join(map(str,incident))
	f.write("Incident particles, " +str(incident)+";\n")
	if type(energy) is list: energy = ','.join(map(str,energy))
	f.write("Mean particle energy, "+str(energy)+";\n")
	f.write("Energy distribution, M, "+str(0)+";\n")
	f.write("Particle fraction, "+str(1.0)+";\n")
	f.write("Dose or fluence, 1, "+str(dose)+";\n")
	f.write("Dose rate, 0.0;\n")
	f.write("Irradiation target, Simple spherical cell model, radius "+str(radius)+" um;\n")
	f.write("Volumes, 0,5,5,5,0,0,0, 1,"+",".join(map(str,geometry[1:4]))+",0,0,0;\n")
	f.write("Chromosome sizes, 46, "+ ", ".join(map(str,[6.1E3/46 for n in range(46)])) +";\n")
	f.write("DNA Density, "+str(DNADensity)+";\n")
	f.write("Cell Cycle Phase, 0;\n")
	f.write("DNA Structure, 0, 1;\n")
	f.write("In vitro / in vivo, 0;\n")
	f.write("Proliferation status, 1,;\n")
	f.write("Microenvironment, 20, 0.01;\n")
	f.write("Damage definition, 0, 0, 10, "+str(bdRange)+", 0, 0;\n")
	f.write("Time, 0;\n")
	f.write("Damage and primary count, "+str(damageCount)+", "+str(primaryCount)+";\n")

	if writeSparse:
		f.write("Data entries, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;\n")
	else:
		f.write("Data entries, "+", ".join(["1"]*14)+";\n")
	f.write("Additional information, Highly simplified DNA damage model for testing. ")
	if writeSparse:
		f.write("Data was generated in minimal output format.;\n")
	else:
		f.write("Data was generated in comprehensive data format.;\n")
	f.write("***EndOfHeader***;\n\n")

# Write hits to file
def writeHits(f, hits):
	for hitSet in hits:
		for event in hitSet:
			f.write('; '.join(map(str,event)))
			f.write(';\n')

# Write dataset to file
def writeToFile(hits, outFile, writeSparse=True, target='Unspecified', geometry=[0,1,1,1],
				DNADensity=-1, bdRange=-1, O2=-1, incident=22, energy=0.1, dose=1.0, 
				function='Gaussian\t0.1\t0.0005', grouping='Single Event'):
	damageCount = sum([len(h) for h in hits])
	primaryCount = sum([sum([min(1,int(event[0][0])) for event in hitSet]) for hitSet in hits])
	with open(outFile,'w') as f:
	 	writeHeader(f, energy, writeSparse, incident, dose, geometry, DNADensity, 
	 				bdRange, damageCount, primaryCount)
	 	writeHits(f, hits)
