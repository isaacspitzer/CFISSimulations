#!/usr/bin/python3

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
import datetime
import glob
import fileinput
import os
#import pyfits
import random
import scipy.stats as stats
import string
import subprocess
import sys
import threading
import time
from UtilityFunctions import ReadVariables

# Get the RA and Dec range for this specific image.
def GetRADecRange(centerRA, centerDec):
	# Calculate the width and height of the image in degrees using the size in pixels and the pixel scale.
	width = float(variables['image_size_x']) * float(variables['pixel_scale']) / 3600.0
	height = float(variables['image_size_y']) * float(variables['pixel_scale']) / 3600.0

	# TODO: Approximating the region of the sky as a flat rectangle...
	minRA = centerRA - width / 2.0
	maxRA = centerRA + width / 2.0
	minDec = centerDec - height / 2.0
	maxDec = centerDec + height / 2.0

	return minRA, maxRA, minDec, maxDec

# This function will generate a new catalog of objects to simulate.  It does not generate any actual images, just the object catalog.
def GenerateCatalog(runName):
    print('GENERATING MOCK CATALOG')
    start = time.time()
    p1 = subprocess.Popen('python ./CreateCatalog.py %s' % runName, shell=True)
    p1.wait()
    end = time.time()
    
    outfile = open('./CatalogGenerationTime.txt', 'w')
    outfile.write('Catalog took %f seconds.' % (end - start))
    outfile.close()

# This function will run in its own thread.
#  It will create a set of simulated images and measure the shapes.
def CreateImages(runName, coreNum, numImagesPerCore, catalogName, ra, dec, headerFile):
    # Generate a seeing and sky noise level for this exposure.
    imageSeeingMean = variables['seeing']
    imageSeeingDispersion = variables['seeing_dispersion']
    imageSeeing = stats.truncnorm(imageSeeingMean - imageSeeingDispersion, imageSeeingMean + imageSeeingDispersion, loc = imageSeeingMean, scale = imageSeeingDispersion).rvs(1)[0]

    skyNoiseLevelMean = variables['sky_' + str(int(variables['sky_brightness']))]
    skyNoiseLevelDispersion = variables['sky_dispersion']
    skyNoiseLevelGauss = stats.truncnorm(skyNoiseLevelMean - skyNoiseLevelDispersion, skyNoiseLevelMean + skyNoiseLevelDispersion, loc = skyNoiseLevelMean, scale = skyNoiseLevelDispersion).rvs(1)[0]
    skyNoiseLevel = float(skyNoiseLevelGauss * variables['gain'] * variables['exposure_time'])

    start = time.time()

    for ccdNum in range(int(variables['numCCDs'])):
        #ccdXOffsetDeg = float(ccd[0]) * float(variables['pixel_scale']) / 3600.0
        #ccdYOffsetDeg = float(ccd[1]) * float(variables['pixel_scale']) / 3600.0
        #minRA, maxRA, minDec, maxDec = GetRADecRange(ra + ccdXOffsetDeg, dec + ccdYOffsetDeg)
	#crpix1 = ccd[0]
	#crpix2 = ccd[1]

        # Run the script that creates the images and wait for it to finish.
        print('python ./CreateImages.py %s %i %i %s %f %f %i %s %f %f' % (runName, coreNum, numImagesPerCore, catalogName, ra, dec, ccdNum, headerFile, imageSeeing, skyNoiseLevel))
        p1 = subprocess.Popen('python ./CreateImages.py %s %i %i %s %f %f %i %s %f %f' % (runName, coreNum, numImagesPerCore, catalogName, ra, dec, ccdNum, headerFile, imageSeeing, skyNoiseLevel), shell=True)
        p1.wait()

        if int(variables['createSextractorCatalog']) == 1:
            CreateCatalog(runName, coreNum, numImagesPerCore, ccdNum)
                
            if int(variables['compareWithCFISExposure']) == 1:
                CompareCatalogs(runName, coreNum, numImagesPerCore, ccdNum)
                CompareWithCFIS(runName, coreNum, numImagesPerCore, ccdNum)

    end = time.time()

    outfile = open('./ImageGenerationTime.txt', 'w')
    outfile.write('Image generation took %f seconds.' % (end - start))
    outfile.close()

# Run Sextractor to create a catalog of detected objects.
def CreateCatalog(runName, coreNum, imageNum, ccdNum):
    print('Running Sextractor on the simulated image.')
    p1 = subprocess.Popen('python ./CreateSextractorCatalogs.py %s %i %i %i' % (runName, coreNum, imageNum, ccdNum), shell=True)
    p1.wait()

# Compare the Sextractor catalog to the known catalog.
def CompareCatalogs(runName, coreNum, imageNum, ccdNum):
    print('Comparing catalogs.')
    p1 = subprocess.Popen('python ./MatchCatalogs.py %s %i %i %i' % (runName, coreNum, imageNum, ccdNum), shell=True)
    p1.wait()
    
# Compare the output of the Sextractor catalog, with the Sextractor catalog of a real CFIS image.
def CompareWithCFIS(runName, coreNum, imageNum, ccdNum):
        p1 = subprocess.Popen('python ./CompareWithCFIS.py %s %i %i %i' % (runName, coreNum, imageNum, ccdNum), shell=True)
        p1.wait()
        
def BulkCompareWithCFIS(runName, coreNum, imageNum):
        p1 = subprocess.Popen('python ./BulkCompareWithCFIS.py %s %i %i' % (runName, coreNum, imageNum), shell=True)
        p1.wait()

# Reads a file containing a list of the RA and Decs of multiple observations.  Returns a list of the [RA, Dec] pairs.
def GetObservations(obsFilename, raIndex, decIndex):
	observations = []

	with open(obsFilename) as infile:
		for line in infile:
			# Skip the first line.
			if '#' in line:
				continue 

			data = line.split()
			headerFile = data[int(variables['headerIndex'])]

			ra = -1.0
			dec = -1.0

			if raIndex > -1 and decIndex > -1:
				c = SkyCoord(data[raIndex] + ' ' + data[decIndex], unit=(u.hourangle, u.deg))
				ra = c.ra.degree
				dec = c.dec.degree

			observations.append([ra, dec, headerFile])

	return observations

def GetObservationsFromHeaders(headerFolder):
	observations = []

	headerFiles = glob.glob(headerFolder + '*.head')

	for headerFile in headerFiles:
		observations.append([-1.0, -1.0, headerFile.split('/')[-1]])

	return observations

# The user is supplying a list of RAs and Decs in the config file.  Parse them.
def GetRADecsFromVariables():
	observations = []
	ras = variables['ras'].split(',')
	decs = variables['decs'].split(',')

	for i in range(len(ras)):
		observations.append([float(ras[i]), float(decs[i])])

	return observations

# Read the CCD layout of the mosaic from a fits header file.
def ReadCCDLayout(ccdLayoutFilename):
	ccds = []
	#with open(ccdLayoutFilename) as infile:
	#	for line in infile:
	#		data = line.split(',')
	#		pixelOffsetX = int(data[0])
	#		pixelOffsetY = int(data[1])
	#		ccds.append([pixelOffsetX, pixelOffsetY])
	hdulist = fits.open(ccdLayoutFilename)
	i = 1
	while i <= variables['numCCDs']:
		crpix1 = float(hdulist[i].header['CRPIX1'])
		crpix2 = float(hdulist[i].header['CRPIX2'])
		ccds.append([crpix1, crpix2])
		i = i + 1

	return ccds

def RewriteVariables(tileName):
	for line in fileinput.input('./Variables.cfg', inplace=1):
		outline = line
		if 'observationName =' in line:
			outline = 'observationName = %s \n' % (tileName)
		elif 'observationsList =' in line:
			outline = 'observationsList = ./ObservingList/%s_ObservingList.txt \n' % tileName
		elif 'dataPath =' in line:
			outline = 'dataPath = ../../Data/Ludo/%s/single_V1.1.0A/ \n' % tileName
		elif 'headerPath =' in line:
			outline = 'headerPath = /home/ispitzer/CFISSimulations/Data/%s/headers_V1.1.0A/ \n' % tileName

		sys.stdout.write(outline)

##############################################################################################
##### BEGINNING OF SCRIPT ####################################################################
##############################################################################################

os.chdir('/home/ispitzer/CFISSimulations/Scripts/')

TILE_NAME = ''
RUN_NUM = 0
if len(sys.argv) > 1:
	TILE_NAME = sys.argv[1]
	RUN_NUM = sys.argv[2]
	RewriteVariables(TILE_NAME)

# Read in the variables from the config file into a dictionary.
variables = ReadVariables()

# Calculate the number of simulations to be created by each core.
numImages = 1
observations = []
print('useObservingList')
print(variables['useObservingList'])

if variables['useObservingList'] == 0:
	numImages = variables['numImages']
	observations = GetRADecsFromVariables()
elif variables['useObservingList'] == 1:
	observations = GetObservations(variables['observationsList'], int(variables['raIndex']), int(variables['decIndex']))
	numImages = len(observations)
else:
	print('Looking to headerPath: %s' % variables['headerPath'])
	observations = GetObservationsFromHeaders(variables['headerPath'])
	print('observations:')
	print(observations)
	numImages = len(observations)

#ccds = [[0, 0]]
#if variables['useCCDLayout'] == 1:
#	ccds = ReadCCDLayout(variables['ccdLayout'])

numCores = int(variables['numCores'])

print(numCores)
print(numImages)
if (numImages < numCores):
	numCores = numImages

#print('numImagesPerCore = int(%i / %i)' % (numImages, numCores))

numImagesPerCore = int(numImages / numCores)

runName = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M") + '_' + str(RUN_NUM)

# If there's already a directory created with this timestamp, wait a minute...
while os.path.isdir('%s%s' % (variables['outdir'], runName)):
    runName = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M") + '_' + str(RUN_NUM)

start = time.time()
catalogName = ''

# Check to see if we want to generate a new catalog or use an existing one.
if int(variables['generateCatalog']) == 1:
	GenerateCatalog(runName)

	#Will need to update this line to match the actual name of the output catalog.
	catalogName = '%s%s/SimulatedCatalog.cat' % (variables['outdir'], runName)
else:
	catalogName = variables['existingCatalog']

threads = []
currentObservation = 0

#for imageNum in range(numImagesPerCore):
print('There are %i observations.' % len(observations))
imageNum = 0

#j = 0
#for i in range(len(observations)):
#	if '2367626' in observations[i][2]:
#		print(observations[i][2])
#		j = j + 1
#
#print(j)
while currentObservation < len(observations):
    # Start all of the threads.
    for coreNum in range(numCores): 
        print('MasterScript Version')
        print(sys.version)   

        if currentObservation >= len(observations):
            continue

        ra, dec, headerFile = observations[currentObservation]
        thr = threading.Thread(target=CreateImages, args=(runName, coreNum, imageNum, catalogName, ra, dec, headerFile))
        threads.append(thr)
        thr.start()
        currentObservation += 1

    # Wait for all of the threads to complete.
    for thr in threads:
        thr.join()

    imageNum += 1

if int(variables['compareWithMultipleCFISExposures']) == 1:
    BulkCompareWithCFIS(runName, coreNum, numImagesPerCore)
elif int(variables['compareWithCFISExposure']) == 1:
    CompareWithCFIS(runName, coreNum, numImagesPerCore)

if len(sys.argv) > 1:
	# Now run MakeCoadds.py to create all of the coadds.
	os.chdir('/cfis/terben/CFIS2000-galsim/CFISCOLLAB_V1.1.0A')
	cfisTile = '%s_r.MP9602.V1.1.0A.swarp.cut.fits' % TILE_NAME
	tileRaDec = TILE_NAME.replace('CFIS_', '')
	inputList = 'inputlist_all_%s.txt' % (tileRaDec)
	cmd = 'python3 MakeCoadds.py %s %s %s %s %s' % (TILE_NAME, cfisTile, RUN_NUM, runName, inputList)
	p1 = subprocess.Popen(cmd, shell=True)
	p1.wait()
	os.chdir('/home/ispitzer/CFISSimulations/Scripts')

end = time.time()
hours = ((end - start) / 60.0) / 60.0
print("The process took %f hours" % hours)
outfile = open('./ProcessDuration.txt', 'w')
outfile.write("The process took %f hours" % hours)
outfile.close()
