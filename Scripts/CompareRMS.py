import datetime
import glob
import math
import matplotlib.pyplot as plt
import numpy
import pyfits
import sys
from astropy.stats import SigmaClip
from photutils import SExtractorBackground
from UtilityFunctions import ReadVariables

def GetHistogram(data, bins):
        print('Making histogram bins')
        print(datetime.datetime.now())
	hist = [0] * len(bins)
	print('Looping... %i bins' % len(bins))
	print(datetime.datetime.now())
	for mag in data:
		ind = numpy.digitize(mag, bins) - 1
		#print('Counts = %i, index = %i' % (int(mag), ind))

		try:
			hist[ind] = hist[ind] + 1
		except:
			print('ind = %i, len(hist) = %i' % (ind, len(hist)))

	return hist

simImageFilename = sys.argv[1]
cfisImageFilename = sys.argv[2]

# Create RMS plots.
simValues = []
#simImageFilenames = glob.glob(simImage)
#for simImageFilename in simImageFilenames:
simImage = pyfits.open(simImageFilename)
simData = simImage[0].data
radiiTotalSim = []
magsTotalSim = []

sigmaClip = SigmaClip(sigma=3.0)
bkg = SExtractorBackground(sigmaClip)
bkgValue = bkg.calc_background(simData)

for row in simData:
	for val in row:
		simValues.append(val - bkgValue)
#../../Data/2074313p.fits
simImage.close()

# Now get the CFIS values.
cfisValues = []
#for cfisCatalogFilename in cfisCatalogFilenames:
cfisImage = pyfits.open(cfisImageFilename)
cfisData = cfisImage[1].data[0:4611, 32:2079]
        
sigmaClip = SigmaClip(sigma=3.0)
bkg = SExtractorBackground(sigmaClip)
bkgValue = bkg.calc_background(cfisData)

for row in cfisData:
	for val in row:
		cfisValues.append(int(val - bkgValue))

cfisImage.close()

print('minSim')
print(datetime.datetime.now())
minSim = min(simValues)
print('maxsim')
print(datetime.datetime.now())
maxSim = max(simValues)
print('minCFIS')
print(datetime.datetime.now())
minCFIS = min(cfisValues)
print('maxCFIS')
print(datetime.datetime.now())
maxCFIS = max(cfisValues)
print(datetime.datetime.now())

# Now make histograms for these pixel values.
print('Getting bins')
print(datetime.datetime.now())
bins = numpy.arange(-150.0, 150.0, 5.0)
print('Making simulation histogram')
print(datetime.datetime.now())
histSim = GetHistogram(simValues, bins)
print('Making CFIS histogram')
print(datetime.datetime.now())
histCFIS = GetHistogram(cfisValues, bins)
print('Plotting')
print(datetime.datetime.now())

print('# Sim values: %i' % len(simValues))
print('# CFIS values: %i' % len(cfisValues))
plt.cla()
plt.clf()
plt.figure(3)
plt.plot(bins, histSim, 'r-', label = 'Simulation', alpha=0.7)
plt.plot(bins, histCFIS, 'b-', label = 'CFIS', alpha=0.7)
plt.ylabel(r'Counts')
plt.xlabel('Pixel Value')
plt.legend(bbox_to_anchor=(0.7, 0.95), loc=2, borderaxespad=0.0)
plt.title('RMS')
#plt.savefig('%sImages/NoiseRMS-%i-%i.png' % (runFolder, coreNum, imageNum))
plt.show()
