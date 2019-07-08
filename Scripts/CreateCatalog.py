import galsim
import math
import numpy as np
import os
import random
from scipy.integrate import quad
from scipy.stats import rv_continuous
import scipy.stats as stats
from shutil import copyfile
import sys
from UtilityFunctions import ReadVariables

# Class for galaxy models.
class Galaxy:
	halfLightRadius = 0.0
	sersicIndex = 0.0
	ellipticity = 0.0
	galType = 0
	bulginess = 0
	magnitude = 0.0
	axisRatio = 0.0

	def __init__(self, hlr, sInd, ell, typ, bul, mag, q):
		self.sersicIndex = sInd
		self.ellipticity = ell
		self.galType = typ
		self.bulginess = bul
		self.magnitude = mag
		self.axisRatio = q
		self.halfLightRadius = hlr

# Get the true shape for an object with the given index.
def GetTruth(trueValues, Ind):
	for truth in trueValues:
		if truth.ind == Ind:
			return (truth.e1, truth.e2)

# Get the number of galaxies with a given magnitude per square degree.
# This function is defined in Fenech Conti et al. 2017.
# The extra factor of 10 comes from the fact that the equation is binned in .1 mag (see the y axis label of Figure 1).
def NumGalaxiesPerSquareDegree(rMag):
	return 10.0 * 10.0 ** (-8.85 + (0.71 * rMag) - (0.008 * rMag * rMag))

# Calculate the number of stars with a given magnitude per square degree.
# This function is based off of Gao, S. (2013).
def NumStarsPerSquareDegree(rMag):
	return 0.1215512948 * math.exp(0.3680992382 * rMag)

# Integrated to correct the image area for the cos(dec) term.
def AreaAdjustment(dec, width):
	return width * math.cos(dec)

# Return the image area in degrees squared.
def GetImageArea(x, y, pixScale):
	areaAS = (x * pixScale) * (y * pixScale)
	areaAM = areaAS / (60.0 * 60.0)
	areaDeg = areaAM / (60.0 * 60.0)
	return areaDeg

# Split up the list of models by magnitude
def PartitionByMagnitude(models, indices):
	idx = 0
	for index in indices:
		sublist = []
		while idx < len(models) and models[idx].magnitude < index:
			sublist.append(models[idx])
			idx = idx + 1

		yield sublist

# Randomly determine an ellipticity using the variables in the config file.
def DrawEllipticity():
	# Draw the ellipticities from a Gaussian with dispersion = variables['ellipticity_dispersion']
	#ellipticity = stats.truncnorm(0.0, upperTrunc, loc = 0.0, scale = variables['ellipticity_dispersion']).rvs(1)

	#return ellipticity
	ellipticity = stats.truncnorm(0.0, upperTrunc, loc = 0.0, scale = variables['ellipticity_dispersion']).rvs(1)
	angle = random.uniform(0.0, 2.0 * np.pi)
	e1 = ellipticity * np.cos(angle)
	e2 = ellipticity * np.sin(angle)

	return e1, e2

# Print a message out to the console, and write it out to the log file as well.
def PrintAndLog(msg, logfile):
	print(msg + '\n')
	Log(msg, logfile)

# Write a message out to the log file without printing it to the screen.
def Log(msg, logfile):
	logfile.write(msg + '\n')

# Read in the catalog of galaxy models.
# This function will return a list of the magnitudes of the models, as well as a list of Galaxy objetcs containing the complete set of model properties.
def ReadModelCatalog(filename):
	mags = []
	models = []
        gemsfile = open(filename, 'r')

        # Skip past the description lines.
        gemsfile.readline()
        line = gemsfile.readline()    
        
        # Keep track of the relevant parameters for all COSMOS galaxies with a reduced Chi^2 of less than 'minChi2'.
        while line:
                data = line.split(',')

                # Uncomment these lines if using COSMOS catalog.	
                #chi2 = float(data[59])
                #cla = int(float(data[21]))

                #if (chi2 < variables['min_chi2'] and cla == variables['galaxy_class']):
                #	hlr = float(data[33]) * variables['catalog_pixel_scale']# * 0.5
                #	sInd = float(data[55])
                #	ell = float(data[43])
                #	typ = int(float(data[64]))
                #	bul = int(float(data[65]))
                #	mag = float(data[6])

                #	if hlr > 0.0 and sInd >= 0.3 and sInd <= 6.0:
                #		mags.append(mag)
                #		models.append(Galaxy(hlr, sInd, ell, typ, bul, mag))
                sInd = 0.0
                hlr = 0.0
                axisRatio = 0.0

                try:
                        sInd = float(data[11])
                        hlr = float(data[9])
                        axisRatio = float(data[13])
                except:
                        line = gemsfile.readline()
                        continue

                if data[3] == 'COSMOS':
                        hlr = hlr * variables['cosmos_pixel_scale'] # COSMOS pixel scale
                else:
                        hlr = hlr * variables['aegis_pixel_scale'] # AEGIS pixel scale

                hlr = hlr
	
                sMag = data[4]
                mag = 0.0

                if sMag.strip() != '':
                        mag = float(sMag)
                else:
                        sMag = -1.0

                starProb = float(data[7])
                objClass = data[5]

		# Filter out stars and other anomolous items.  GalSim can only handle Sersic indices between 0.3 and 6.2 inclusive (galsim-developers.github.io/GalSim/classgalsim_1_1base_1_1_sersic.html)
                if mag > 0.0 and hlr > 0.0 and hlr < 3.0 and sInd >= 0.3 and sInd <= 6.2 and starProb < 0.2 and objClass.find('Star') == -1 and objClass.find('QSO') == -1:
                        mags.append(mag)
                        models.append(Galaxy(hlr, sInd, 0.0, 0, 0, mag, axisRatio))			

                line = gemsfile.readline()

        return (mags, models)


def AddGalaxyToCatalog(mag, hlr, sersicIndex, e1, e2, galCount, numRotations, catalogFile, catWidth, catHeight):
	currentGalaxyRotation = 0

	# Pick an initial angle for this galaxy.
	#angle = random.uniform(0.0, 2.0 * np.pi)
	angleToRotate = variables['rotation_angle']
	angleToRotate = angleToRotate * (np.pi / 180.0)
	shape = np.matrix('%f %f' % (e1, e2))

	# Add the correct number of rotations for this galaxy.
	# Generate new positions and orientations for this galaxy.
	for j in range(numRotations):
		posX = random.uniform(0, catWidth) + variables['minRA']
                #posY = random.uniform(0, catHeight) + variables['minDec']
		# The declination should be drawn from a cosine distribution so that the regions
		#  closer to the poles don't get overcrowded.
		posY = math.degrees(cosDist.rvs(size=1)[0])

		# Rotate the e1, e2 values.
		cos2Phi = np.cos(2.0 * angleToRotate * j)
		sin2Phi = np.sin(2.0 * angleToRotate * j)
		rotation = np.matrix('%f %f;%f %f' % (-1.0 * cos2Phi, -1.0 * sin2Phi, sin2Phi, -1.0 * cos2Phi))
		newShape = shape.dot(rotation)
		
		newE1 = newShape.item(0)
		newE2 = newShape.item(1)

		# 'Gal#, X, Y, r-Mag, HLR,  Sersic Index, e1 (intrinsic), e2 (intrinsic) \n'
		catalogFile.write('%i %f %f %f %f %f %f %f \n' % (galCount + currentGalaxyRotation, posX, posY, mag, hlr, sersicIndex, newE1, newE2))

		#angle += angleToRotate
		currentGalaxyRotation += 1

def AddStarToCatalog(starNum, mag, e1, e2, catalogFile, catWidth, catHeight):
	# TODO: If we want to implement picking seeing from a Gaussian distribution, this line will need to be updated.
	imageSeeing = variables['seeing']

        posX = random.uniform(0, catWidth) + variables['minRA']
        posY = math.degrees(cosDist.rvs(size=1)[0])

        catalogFile.write('%i %f %f %f %f %f %f %f\n' % (starNum, posX, posY, mag, imageSeeing, -1, e1, e2))

def cosFunction(dec):
	return (np.radians(variables['maxRA']) - np.radians(variables['minRA'])) * np.cos(dec)

def CalculateCatalogArea():
	return quad(cosFunction, np.radians(variables['minDec']), np.radians(variables['maxDec']))[0]

###############################################################################################################
###     BEGINNING OF SCRIPT     ###############################################################################
###############################################################################################################

runName = sys.argv[1]

# Create the folder for the images.
outfolder = '../Runs/%s/Images/' % (runName)
if not os.path.isdir(outfolder):
	try:
		os.makedirs(outfolder)
	except:
		print('Output folder was created by another process.')

variables = ReadVariables()

#####################################################################
# The following code block allows me to draw random number from a cos
#  distribution.  That way, high dec regions to not get overcrowded
#  due to the Cos(dec) shringking of the regions near the poles.
#####################################################################

# This function is integrated with the dec bounds to get the normalization constant.
def UnnormalizedPDF(x):
	return np.cos(x)

limLower = math.radians(variables['minDec'])
limUpper = math.radians(variables['maxDec'])
norm = quad(UnnormalizedPDF, limLower, limUpper)[0]

# Draw random values from a bounded cosine distribution.
class CosDistribution(rv_continuous):
	def _pdf(self, x):
		return np.cos(x) / norm

#####################################################################
# End of cos distribution code block.
#####################################################################

catalogFilename = '../Runs/%s/SimulatedCatalog.cat' % runName

logfile = open('../Runs/%s/log-CreateCatalog.txt' % (runName), 'w')
msg = 'Creating output folders'
PrintAndLog(msg, logfile)

# Make a local copy of the config file for future reference.
if not os.path.isfile('../Runs/%s/Variables.cfg' % runName):
	copyfile('Variables.cfg', '../Runs/%s/Variables.cfg' % runName)

msg = 'Setting up Variables'
PrintAndLog(msg, logfile)

cosDist = CosDistribution(a = limLower, b = limUpper)

lowerTrunc = (-1.0 * variables['shape_cutoff'] - variables['shape_mean']) / variables['ellipticity_dispersion']
upperTrunc = (variables['shape_cutoff'] - variables['shape_mean']) / variables['ellipticity_dispersion']

# Since COSMOS data was acquired with HST, I need to convert the counts so that they make sense for CFIS data.
fluxScaling = (variables['telescope_diameter']**2 / 2.4**2 * (1.-0.33**2))

imageAreaDegrees = CalculateCatalogArea() * (180.0 / np.pi) * (180.0 / np.pi) #(variables['maxRA'] - variables['minRA']) * (variables['maxDec'] - variables['minDec'])

# quad(NumGalaxiesPerSquareDegree, 0, variables['mag_max'])[0]
#imageWidth = (math.radians(variables['maxRA']) - math.radians(variables['minRA']))
#imageArea = quad(AreaAdjustment, math.radians(variables['minDec']), math.radians(variables['maxDec']), args=(imageWidth))[0]
#PrintAndLog('Adjusted catalog area (sq rad): %f' % imageArea, logfile)
#imageAreaDegrees = imageArea * (180.0 / np.pi) * (180.0 / np.pi)

#PrintAndLog('Adjusted catalog area (sq deg): %f' % imageAreaDegrees, logfile)

catWidth = (variables['maxRA'] - variables['minRA'])
catHeight = (variables['maxDec'] - variables['minDec'])

msg = 'Reading COSMOS catalog.'
PrintAndLog(msg, logfile)

# Read in parameters from COSMOS catalog.
mags, models = ReadModelCatalog('../COSMOS/asu.csv')

msg = 'Done reading COSMOS catalog. min mag: %f, max mag: %f, # Models: %f' % (min(mags), max(mags), len(mags))
PrintAndLog(msg, logfile)

models.sort(key = lambda x: x.magnitude, reverse = False)

partitionMags = np.arange(11.0, 30.0, variables['mag_step'])
modelMagLists = list(PartitionByMagnitude(models, partitionMags))

leftovers = 0
galCount = 0
numRotations = int(variables['num_rotations'])
galDensityBoostFactor = 1.0
if variables['gal_density_target'] > 0.0:
	galDensityBoostFactor = variables['gal_density_target'] / (quad(NumGalaxiesPerSquareDegree, 0, variables['mag_max'])[0] / 3600.0) #This constant is the integral of the number density function from 0->mag_max (Fenech Conti)

catalogFile = open(catalogFilename, 'w')
catalogFile.write('Gal#, X, Y, r-Mag, HLR,  Sersic Index, e1 (intrinsic), e2 (intrinsic) \n')

loopCount = 0

# Loop through each magnitude bin, adding the appropriate number of galaxies given the size of the field.
for i in range(len(modelMagLists)):
	loopCount += 1

	# If there are no galaxy models in this magnitude range, add the number we would have added to the catalog to 'leftovers', and continue.
	if len(modelMagLists[i]) == 0:
		currentMag = partitionMags[loopCount - 1] + (variables['mag_step'] / 2.0)
		numGalaxies = NumGalaxiesPerSquareDegree(currentMag) * imageAreaDegrees * galDensityBoostFactor * variables['mag_step']
		msg = 'No models in range %f - %f, currentMag = %f' % (partitionMags[i], partitionMags[i] + variables['mag_step'], currentMag)
		PrintAndLog(msg, logfile)
		msg = 'NumGalaxiesPerSquareDegree = %f, imageArea = %f, galDensityBoostFactor = %f, mag_step = %f' % (NumGalaxiesPerSquareDegree(currentMag), imageAreaDegrees, galDensityBoostFactor, variables['mag_step'])
		PrintAndLog(msg, logfile)
		msg = 'NumGalaxies = %f, leftovers = %f' % (numGalaxies, leftovers)
		PrintAndLog(msg, logfile)
		leftovers = leftovers + numGalaxies
		continue
	# If the magnitude range is above the limit, skip them.
	if modelMagLists[i][0].magnitude > variables['mag_max']:
		msg = 'Skipping bin with mag %f.' % modelMagLists[i][0].magnitude
		PrintAndLog(msg, logfile)
		continue

	# Calculate the appropriate number of galaxies for this magnitude range, get an appropriate set of parameters for them, and add them to the catalog.
	# Calculate the number of galaxies of this magnitude to add to the image.
	numGalaxies = NumGalaxiesPerSquareDegree((modelMagLists[i][0].magnitude + modelMagLists[i][-1].magnitude) / 2.0) * imageAreaDegrees * galDensityBoostFactor * variables['mag_step']

	msg = 'Creating %i galaxies with magnitudes between %f and %f + %f leftovers.' % (numGalaxies, modelMagLists[i][0].magnitude, modelMagLists[i][-1].magnitude, leftovers)
	PrintAndLog(msg, logfile)

	# If the previous magnitude bin had issues, add the number of galaxies from that bin to this one to compensate.
	numGalaxies = numGalaxies + leftovers

	# If we want multiple rotations of the same galaxy, reduce the number of models to draw.
	numUniqueGalaxies = int(numGalaxies / variables['num_rotations']) 
	print('Number of unique galaxy models to use: %i' % numUniqueGalaxies)

	# If there were 5 galaxies to add, but we have 4 rotations, then 1 galaxy is not added.  Add it to the next bin.  This is really only important for low magnitude galaxies where the counts are low.
	leftovers = numGalaxies % int(variables['num_rotations'])

	# Draw properties for the galaxies and add them to the catalog.
	for galNum in range(numUniqueGalaxies):
		# Randomly pick a model from the appropriate magnitude range.
		galInd = random.uniform(0, len(modelMagLists[i]) - 1)
		model = modelMagLists[i][int(galInd)]

		# Draw the ellipticities from a Gaussian with dispersion = variables['ellipticity_dispersion']		
		e1, e2 = DrawEllipticity()

		# Because half-light radius depends on axis ratio, we need to figure out the axis ratio first to correct the HLR.
		#angle = 0.0
		#e1 = ellip * np.cos(angle)
		#e2 = ellip * np.sin(angle)
		axisRatio = 1.0 / math.exp(galsim.Shear(g1 = e1, g2 = e2).eta)
		hlr = model.halfLightRadius * math.sqrt(axisRatio)

		AddGalaxyToCatalog(model.magnitude, hlr, model.sersicIndex, e1, e2, galCount, numRotations, catalogFile, catWidth, catHeight)
		galCount += numRotations

# Now add the stars.

msg = 'Adding stars to the catalog.'
PrintAndLog(msg, logfile)

starCount = 0
leftovers = 0.0

starDensityBoostFactor = 1.0
if variables['star_density_target'] > 0.0:
	starDensityBoostFactor = variables['star_density_target'] / (quad(NumStarsPerSquareDegree, 0, variables['mag_max'])[0] / 3600.0)

for mag in partitionMags:

	# If we're past the requested magnitude, stop making stars.
	if mag > variables['mag_max']:
		break
            
        numStars = NumStarsPerSquareDegree(mag) * imageAreaDegrees * starDensityBoostFactor * variables['mag_step']
        remainder = numStars % 1.0
        leftovers = leftovers + remainder
        numStars = int(numStars + leftovers)

        if leftovers > 1.0:
                leftovers = leftovers % 1.0

        msg = 'Creating %i stars with magnitude %f.' % (numStars, mag)
        Log(msg, logfile)

        for starNum in range(numStars):
                AddStarToCatalog(starCount, mag, 0.0, 0.0, catalogFile, catWidth, catHeight)
		starCount += 1

catalogFile.close()
