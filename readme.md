
# Medras-MC

This is a Monte Carlo implementation of the MEDRAS (Mechanistic DNA Repair and Survival) model. This model simulates the response of cells to DNA damage following ionising radiation exposure, including DNA repair kinetics, accuracy of repair, and resulting cell fate. 

Descriptions of the model and an analytical implementation have been published in (McMahon 2016) and (McMahon 2017). This code presents a Monte Carlo based implementation of the model, which enables the accurate simulation of arbitrary DNA damage distributions, and a more detailed assessment of the resulting DNA alterations which follow from irradiation.

This code is designed to read in DNA damage data stored in the Standard for DNA Damage (Schuemann 2019), which enables the ready exchange of DNA damage data between different simulation codes. 

This distribution consists of two components - damagegenerator, which can be used to generate DNA damage distributions based on a simple energy-deposition based model, and repairanalysis, which analyzes repair in input SDD files. These are described briefly below.

## damagegenerator


This code implements a simplified model of DNA damage following ionising radiation exposure. This assumes that the average yield of damage within a given volume is proportional to the average energy deposited within that volume, and that these events are randomly distributed. For X-rays this is modelled as simple random damage throughout the nucleus, while for ions this is modelled based on the energy distribution around tracks passing through the nucleus. Damage is flagged as either 'simple' or 'complex' based on a simple probability.

`damagegenerator` can be used to create a damage distribution as follows:

```py
# Add damagegenerator to path as preferred
from damagegenerator import damageModel

damageModel.basicXandIon()  # Generate an illustrative dataset
damageModel.generateExposure(10.0, 4.58, 1, 1, 10) # Generate event for 10 MeV proton
```

The full argument list for generateExposure is:

- Energy (in MeV)
- LET (in keV/μm)
- Dose (in Gy)
- Particle atomic number (Z)
- Number of exposures to simulate
- Nucleus radius (μm, default 4.32)
- Chromosome count (46)
- Additional Info (string)

These commands will generate SDD files containing the resulting DSB distribution. By default, only the spatial information is meaningful in these data, so the files by default are output in a minimal SDD format, excluding fields which relate to data which are not simulated in this model. Optionally, 'writeSparse' in damageModel.py can be set to false, which will write a fully-detailed SDD file with synthetic data for other features. This may be useful for testing purposes, but should not be taken as biologically meaningful.

## repairanalysis

This code implements the MEDRAS repair model for a set of input SDD files. In this approach, the code extracts DSB distributions per exposure from the SDD files, and probabalistically simulates their repair, modelling the stochastic rejoining of different DSB ends. This simulation tracks the yield of misrepaired breaks, and is capable of further sub-analysis including yields of chromosome aberrations, fragment loss, and spatial distributions of misrepair. Several example analyses are available by default. 

`repairanalysis` can be called as follows:

```py
# Add repairanalysis to path as preferred
from repairanalysis import medrasrepair

sddpath = "/path/to/sdd/folder/"
print('Fidelity analysis')
medrasrepair.repairSimulation(sddpath,'Fidelity')

print('\n\nSeparation analysis')
medrasrepair.repairSimulation(sddpath,'Separation')

print('\n\nMisrepair spectrum analysis')
medrasrepair.repairSimulation(sddpath,'Spectrum')
```

Each of the three options presents a different sub-analysis of the results of repair, summarised below. Options to modify the output of each of these processes are available at the beginning of medrasrepair.py. Common parameters allow for the number of exposures per file to be altered, as well as the number of times the simulation should be repeated for each exposure. 

### Fidelity

Fidelity (the default approach) summarises the number of breaks, probability of misrepair, and rate of inter-chromosome events when repair is simulated. These will be listed for each individual row as the model proceeds, and then summary data will be presented showing average rates across the whole dataset, as well as data on the kinetics of misrepair. This can either be the time taken to clear DNA DSB repair foci (if `addFociDelay` is true in medrasrepair.py) or simply the time for physical break rejoining (if `addFociDelay` is false). The former option provides better agreement with foci data such as γH2AX or 53BP1 analysis, while the latter option provides better agreement with physical measures of DSB repair such as Pulsed Field Gel Electrophoresis. 

### Separation

Separation summarises the spatial separation between DSBs which misrepair. This is presented as the probability density of misrepair events occurring at different separations through the entire input SDD file, across the entire nucleus volume. This can be useful for identifying regions which provide the greatest contribution to misrepair events.

### Spectrum

The misrepair spectrum analysis calculates the types of misrepair event which result from a given exposure. For the requested timepoint (which can be altered by modifying `simulationLimit` in medrasrepair.py), it calculates the number of residual breaks, the number of misrepaired breaks, 'large' misrepairs (involving movement of more than 3 mega-basepair (MBP) by default), and the number of inter-chromosome misrepairs. It can then summarise the number of chromosome aberrations in several ways, depending on user options. 

If `analyzeAberrations.printChromosomeTypes` is True (default False), it will produce a full summary of the post-repair chromosomes, including the number of chromosomes with 1 or more than 1 chromosome junctions (i.e. inter-chromosome misrepairs), and the number of acentric, normal and multi-centric chromosomes with different structures. 

If `analyzeAberrations.printAberrations` is True (default True), it will produce a simplified summary of chromosome aberrations, based on standard counting methods - dicentrics and multi-centrics; rings; excess linear fragments; and total aberrations. This has been designed to approximately correspond to traditional solid stain (e.g. Giemsa) counting of aberrations.

If `analyzeAberrations.printViable` is True (default True), it will also report if a cell is expected to be viable or not - that is, is it free from aberrations which would prevent replication. Note that this does not correspond directly to survival in all cells, as some cells have active arrest or apoptotic pathways which can lead to loss of viability based on overall levels of damage, independent of the eventual misrepair.

Finally, if  `analyzeAberrations.printDNALoss` is True (default False), the amount of DNA in small fragments initially and at the final timepoint are presented (in MBP), which may be compared with e.g. quantities of DNA released in PFGE studies.

In addition, if `doPlot` is set to true, the model will also generate an mFISH-style plot of the resulting chromosome distribution for each simulated repair, to give a visual illustration of the resulting chromosome alterations.

## Requirements

This code is written in python3, and requires the following libraries:

- numpy
- scipy
- openpyxl
- matplotlib

## Contacts

For questions/comments/bug reports, please contact stephen.mcmahon (at) qub.ac.uk

## References

[McMahon2017]	McMahon, S. J., McNamara, A. L., Schuemann, J., Paganetti, H., & Prise, K. M. (2017). A general mechanistic model enables predictions of the biological effectiveness of different qualities of radiation. Scientific Reports, 7(1), 688. http://doi.org/10.1038/s41598-017-10820-1

[McMahon2016]	McMahon, S. J., Schuemann, J., Paganetti, H., & Prise, K. M. (2016). Mechanistic Modelling of DNA Repair and Cellular Survival Following Radiation-Induced DNA Damage. Scientific Reports, 6, 33290. http://doi.org/10.1038/srep33290

[McMahon2021]   McMahon, S. J. & Prise, K. M. (2021). A Mechanistic DNA Repair and Survival Model (Medras): Applications to Intrinsic Radiosensitivity, Relative Biological Effectiveness and Dose-Rate. Frontiers in Oncology, 11, 689112. https://doi.org/10.3389/fonc.2021.689112

[Schuemann2019]	Schuemann J., McNamara A. L., Warmenhoven J., et al. (55 authors). (2019). A new Standard DNA Damage (SDD) data format. Radiation Research; 191(1):76. PMID: 30407901. PMCID: PMC6407706. https://doi.org/10.1667/RR15209.1
