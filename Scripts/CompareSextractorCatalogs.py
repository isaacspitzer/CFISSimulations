from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import numpy as np
import sys
from UtilityFunctions import ReadVariables

def PartitionByMagnitude(models, indices):
	idx = 0
	for index in indices:
		sublist = []
		while idx < len(models) and models[idx] < index:
			sublist.append(models[idx])
			idx = idx + 1
		yield sublist

simulationCatalog = sys.argv[1]
cfisCatalog = sys.argv[2]
#imageNum = int(sys.argv[3])

variables = ReadVariables()

pixelScale = variables['pixel_scale'] #0.187

#runFolder = '../Runs/%s/' % runName
#galsimCatalog = runFolder + 'Catalog-%i-%i.txt' % (coreNum, imageNum)
#sextractorCatalog = runFolder + 'FullImage-%i-%i-0Cat.txt' % (coreNum, imageNum)
#imageName = runFolder + 'FullImage-%i-%i-0.fits' % (coreNum, imageNum)

gsCat = open(simulationCatalog, 'r')

x1s = []
x2s = []
y1s = []
y2s = []
mags = []

for line in gsCat:
	if '#' in line:
		continue

	data = line.split()
	#x = (float(data[int(variables['sextractorXPosIndex'])]) * pixelScale) / 3600.0
	#y = (float(data[int(variables['sextractorYPosIndex'])]) * pixelScale) / 3600.0

	#x1s.append(x)
	#y1s.append(y)
	flux = float(data[int(variables['sextractorFluxIndex'])])
	zpMag = (-2.5 * np.log10(flux / float(variables['exposure_time']))) + (float(variables['zero_point_r']) * float(variables['gain']))

	mags.append(zpMag)

sCat = open(cfisCatalog, 'r')
numLines = sum(1 for line in open(cfisCatalog))
lineLimit = numLines / 40
lineCount = 0
sMags = []
#sMags2 = []

for line in sCat:
	lineCount = lineCount + 1

	if line.startswith('#') or lineCount > lineLimit:
		continue
	
	data = line.split()
	#x = (float(data[int(variables['sextractorXPosIndex'])]) * pixelScale) / 3600.0
	#y = (float(data[int(variables['sextractorYPosIndex'])]) * pixelScale) / 3600.0
	#counts = float(data[0])
	#zpMag = -2.5 * np.log10(counts / float(variables['exposure_time'])) + float(variables['zero_point_r'])
	#sexMag = float(data[int(variables['sextractorMagIndex'])])
	flux = float(data[int(variables['sextractorFluxIndex'])])
	zpMag = -2.5 * np.log10(flux / float(variables['exposure_time'])) + float(variables['zero_point_r'])

	#x2s.append(x)
	#y2s.append(y)
	#sMags.append(sexMag)
	sMags.append(zpMag)

#c1 = SkyCoord(ra = np.asarray(x1s) * u.degree, dec = np.asarray(y1s) * u.degree)
#c2 = SkyCoord(ra = np.asarray(x2s) * u.degree, dec = np.asarray(y2s) * u.degree)

#idx, d2d, d3d = c2.match_to_catalog_sky(c1)

#detectedMags = []
#sexMags = []

#i = 0

#for ind in idx:
#	detectedMags.append(mags[ind])
#	sexMags.append(sMags[i])

#	i = i + 1

#plt.figure(0)
#plt.plot(detectedMags, sexMags, 'r,')
#plt.xlabel('Actual Magnitude')
#plt.ylabel('Sextractor Magnitude')
#plt.xlim([18.0, float(variables['mag_max']) + 1.0])
#plt.ylim([18.0, float(variables['mag_max']) + 1.0])
#plt.savefig('%sImages/MagnitudeComparison-%i-%i.png' % (runFolder, coreNum, imageNum))
#plt.close(0)

#idx = np.sort(idx)
#detectedMags = []
#for ind in idx:
#	detectedMags.append(mags[ind])

simDetectedCount = 0
cfisDetectedCount = 0

for mag in mags:
	#if mag <= 25.0:
	simDetectedCount = simDetectedCount + 1

for mag in sMags:
	#if mag <= 25.0:
	cfisDetectedCount = cfisDetectedCount + 1

imageAreaAM = ((float(variables['image_size_x']) * pixelScale) * (float(variables['image_size_y']) * pixelScale)) / 3600.0
simSourceDensity = simDetectedCount / imageAreaAM
cfisSourceDensity = cfisDetectedCount / imageAreaAM
print "The density of r < 25 galaxies in the simulation is %f." % simSourceDensity
print "The density of r < 25 galaxies in the CFIS field is %f." % cfisSourceDensity

#partitionMags = np.arange(1.0, 40.0, 0.1)
#simPartitions = list(PartitionByMagnitude(mags, partitionMags))
#cfisPartitions = list(PartitionByMagnitude(sMags, partitionMags))

#print(mags)
#print(sMags)

#xs = []
#y1s = []
#y2s = []
#i = 0

#for partitionMag in partitionMags:
#	xs.append(partitionMag)
#	y1s.append(len(simPartitions[i]))
#	y2s.append(len(cfisPartitions[i]))
#	i = i + 1

#histSim, binEdgesSim = np.histogram(mags, 30)
#histCFIS, binEdgesCFIS = np.histogram(sMags, 30)

plt.cla()
plt.clf()
plt.figure(1)
#plt.semilogy(xs, y1s, 'r-', label='Simulated Galaxies')
#plt.semilogy(xs, y2s, 'b--', label='CFIS Galaxies')
#plt.plot(xs, y1s, 'r-', label='Simulated Galaxies')
#plt.plot(xs, y2s, 'b--', label='CFIS Galaxies')
plt.hist(mags, facecolor='red', alpha=0.5, label='Simulated Galaxies')
plt.hist(sMags, facecolor='blue', alpha=0.5, label='CFIS Galaxies')
plt.xlabel('Magnitude')
plt.ylabel('# of Galaxies')
plt.legend(bbox_to_anchor=(0.45, 1))
#plt.xlim([18.0, float()])
#plt.savefig('%sImages/SextractorCompleteness-%i-%i.png' % (runFolder, coreNum, imageNum))
plt.show()
plt.close(1)
