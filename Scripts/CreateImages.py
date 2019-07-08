from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
import galsim
import math
import numpy
import os
import random
from shutil import copyfile
import sys
import glob
import pyfits
import scipy.stats as stats
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

# Return the image area in degrees squared.
def GetImageArea(x, y, pixScale):
	areaAS = (x * pixScale) * (y * pixScale)
	areaAM = areaAS / (60.0 * 60.0)
	areaDeg = areaAM / (60.0 * 60.0)
	return areaDeg

# Print a message out to the console, and write it out to the log file as well.
def PrintAndLog(msg, logfile):
	print(msg + '\n')
	Log(msg, logfile)

# Write a message out to the log file without printing it to the screen.
def Log(msg, logfile):
	logfile.write(msg + '\n')
    
# Converts the 'e1s' and 'e2s' variables into lists of paired shear values.
def GenerateShearValues(e1s, e2s):
        e1sStr = e1s.split(',')
        e2sStr = e2s.split(',')
        
        shears = []
        
        for i in range(len(e1sStr)):
                shears.append([float(e1sStr[i]), float(e2sStr[i])])
            
        return shears

# Takes a Galsim PositionD in radians and returns the RA and Dec in degrees.
def PositionRadiansToDegrees(pos):
	ra = (pos.ra / galsim.radians) * (180.0 / numpy.pi)
	dec = (pos.dec / galsim.radians) * (180.0 / numpy.pi)
	return ra, dec

# Takes a set of model parameters, generates a galaxy and adds it to the provided image.
def AddGalaxyToImage(ra, minRA, dec, minDec, mag, hlr, sInd, e1, e2, galCount, psf, fullImages, psfImage, truthFile, note, bcoImage, wcs):
        galaxy = None

	# Calculate the location of the galaxy in pixel space.
	worldPosition = galsim.CelestialCoord(ra=ra * galsim.degrees, dec=dec * galsim.degrees)
	
	try:
		imagePosition = wcs.toImage(worldPosition)
	except:
		PrintAndLog("RA: %f, Dec %f" % (ra, dec), logfile)
		return

	#if imagePosition.x < 0.0 or imagePosition.x > variables['image_size_x'] or imagePosition.y < 0.0 or imagePosition.y > variables['image_size_y']:
	#	return None
	#xPosDeg = ra - minRA
	#yPosDeg = dec - minDec
	#xPosPix = xPosDeg * 3600.0 / variables['pixel_scale']
	#yPosPix = yPosDeg * 3600.0 / variables['pixel_scale']

        if note == 'BCO':
                galaxy = galsim.Sersic(n = sInd, scale_radius = hlr)
	elif note == 'Star':
		#galaxy = psf
		galaxy = galsim.Moffat(beta = 3, fwhm = imageSeeing)
        else:
                galaxy = galsim.Sersic(n = sInd, half_light_radius = hlr)
		#galaxy = wcs.toImage(galaxy, image_pos = imagePosition)

        # Calculate the number of counts by using the zero point.  m = -2.5*log_10(counts) + zero_point
        # Zero point is in e- per sec; multiply by gain to convert to ADU per sec.
        counts = variables['exposure_time'] * 10.0 ** ((mag - (variables['zero_point_r'] * variables['gain'])) / -2.5)
        galaxy = galaxy.withFlux(counts)
        
	# Give the galaxy it's intrinsic shape.
        gal_shape = galsim.Shear(g1=e1, g2=e2)
	galaxy = galaxy.shear(gal_shape)

        for i in range(numCopiesOfImage):
		final = None

                # Apply the shear (the actual shear, not the intrinsic shape) and convolve with the PSF.
		if note != 'Star':
	                shearedGalaxy = galaxy.shear(g1 = shears[i][0], g2 = shears[i][1])
			final = galsim.Convolve([psf, shearedGalaxy], gsparams = big_fft_params)
		else:
			# If the object is a star, it doesn't need to be sheared or convolved.
			final = galaxy
			
                stamp = None
                # Draw the convolved galaxy onto the image.			
                try:
			if note != 'Star':
	                        stamp = final.drawImage(scale = variables['pixel_scale'])
			else:
			        stamp = final.drawImage(scale = variables['pixel_scale'], nx=128, ny=128)
                except:
                        errFile = open('./BadModels.txt', 'a')
                        msg = 'Memory error. Magnitude: %f, HLR: %f, Sersic Index: %f' % (mag, hlr, sInd)
                        PrintAndLog(msg, logfile)
                        errFile.write(msg)
                        errFile.close()
                        break

		# Place the object on the image.
                #stamp.setCenter(xPosPix, yPosPix)
		stamp.setCenter(imagePosition.x, imagePosition.y)
		
                bounds = stamp.bounds & fullImages[i].bounds

		try:
	                fullImages[i][bounds] += stamp[bounds]
		except:
			# If the galaxy stamp is too far outside the frame, this will fail.  Just return None.
			PrintAndLog('Galaxy too far out of frame.  Not adding.', logfile)
			return None

                if i == 0:
                        # Create the image and draw the PSF.  Only need to do this once for each set of shear variations.
                        psfStamp = galsim.ImageF(stamp.bounds)
                        try:
                                psf.drawImage(psfStamp, wcs = wcs)

				if int(variables['write_psf_file']) == 1:
                                	psfImage[bounds] += psfStamp[bounds]
                        except:
                                msg = "Missing a PSF"
                                PrintAndLog(msg, logfile)
                                break
                    
                        if note == 'BCO':
                                bounds = stamp.bounds & bcoImage.bounds
				if int(variables['write_bco_file']) == 1:
	                                bcoImage[bounds] += stamp[bounds]
			
                        # Write the galaxy to the output catalog.  Since multiple images zill have the same set of catalogs, only write out once.
                        truthFile.write('%i %f %f %f %f %f %f %f %s %f %f\n' % (galCount, imagePosition.x, imagePosition.y, mag, hlr, sInd, e1, e2, note, ra, dec))

# Get the core # and number of images per core from the command line arguments.
runName = sys.argv[1]
coreNum = int(sys.argv[2])
imageNum = int(sys.argv[3])
catalogFilename = sys.argv[4]
ra = float(sys.argv[5])
dec = float(sys.argv[6])
ccdNum = int(sys.argv[7])
headerFile = sys.argv[8]
imageSeeing = float(sys.argv[9])
skyNoiseLevel = float(sys.argv[10])

# We don't want to just draw galaxies that are 'centered' in the frame.  We also want galaxies that are centered outside the frame, but which overlap with the borders of the frame and bleed into the image.  This setting dictates how far outisde of the frame (in degrees) a galaxy can be and still be added to the image.
OVERLAP_CHECK_DEG = 0.0015

# Specify an output folder for all of the generated images.
coreName = coreNum
outfolder = '../Runs/%s/Images/' % (runName)

# Create the folder for the images.
if not os.path.isdir(outfolder):
	try:
		os.makedirs(outfolder)
	except:
		print('Output folder was created by another process.')

logfile = open('../Runs/%s/log-CreateImages-%i-%i.txt' % (runName, imageNum, ccdNum), 'w')
msg = 'Creating output folders'
PrintAndLog(msg, logfile)

variables = ReadVariables()

# Make a local copy of the config file for future reference.
if not os.path.isfile('../Runs/%s/Variables.cfg' % runName):
	copyfile('Variables.cfg', '../Runs/%s/Variables.cfg' % runName)

msg = 'Setting up Variables'
PrintAndLog(msg, logfile)

big_fft_params = galsim.GSParams(maximum_fft_size = int(variables['maximum_fft_size']))

if int(variables['random_seed']) >= 0:
        rng = galsim.BaseDeviate(int(variables['random_seed']) + 1)
        ud = galsim.UniformDeviate(int(variables['random_seed']) + 1)
else:
        rng = galsim.BaseDeviate(random.randint(0, 5000000))
        ud = galsim.UniformDeviate(random.randint(0, 5000000))

# Since COSMOS data was acquired with HST, I need to convert the counts so that they make sense for CFIS data.
fluxScaling = (variables['telescope_diameter']**2 / 2.4**2 * (1.-0.33**2))

models = []
mags = []

msg = 'Reading input catalog.'
PrintAndLog(msg, logfile)

# Define our simple world coordinate system.
cfisImage = None
if variables['useObservingList'] == 1:
	cfisImage = fits.open(variables['dataPath'] + headerFile)
else:
	cfisImage = fits.open(variables['wcs_image'])

cfisHeader = cfisImage[ccdNum+1].header

# Open the header and rewrite the telescope pointing.
gsHeader = galsim.FitsHeader(header = cfisHeader)
# TODO: Approximating as flat rectangle.
c = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
cHMS = c.to_string('hmsdms')
raString = cHMS.split()[0].replace('h', ':').replace('m', ':').replace('s', '')
decString = cHMS.split()[1].replace('d', ':').replace('m', ':').replace('s', '').replace('+', '')
gsHeader['ra'] = raString
gsHeader['dec'] = decString
gsHeader['crval1'] = ra
gsHeader['crval2'] = dec
gsHeader['ra_deg'] = ra
gsHeader['dec_deg'] = dec
gsHeader['datasec'] = '[%i:%i,1:%i]' % (variables['calibration_columns'] + 1, variables['image_size_x'] - variables['calibration_columns'], variables['image_size_y'] - variables['calibration_columns'])
gsHeader['fscale'] = variables['fscale']

wcs = galsim.wcs.readFromFitsHeader(gsHeader)[0]
PrintAndLog("WCS is a %s" % wcs.__class__.__name__, logfile)

#  Want to find RA, Dec limits of the image.  Get the RA/Decs of each corner and find the min and max.
pos1 = wcs.toWorld(galsim.PositionD(1.0, 1.0))
pos2 = wcs.toWorld(galsim.PositionD(1.0, variables['image_size_y']))
pos3 = wcs.toWorld(galsim.PositionD(variables['image_size_x'], variables['image_size_y']))
pos4 = wcs.toWorld(galsim.PositionD(variables['image_size_x'], 1.0))

# Convert the radians to degrees.
ra1, dec1 = PositionRadiansToDegrees(pos1)
#ra1 = ra1 * math.cos(math.radians(dec1))
ra2, dec2 = PositionRadiansToDegrees(pos2)
#ra2 = ra2 * math.cos(math.radians(dec2))
ra3, dec3 = PositionRadiansToDegrees(pos3)
#ra3 = ra3 * math.cos(math.radians(dec3))
ra4, dec4 = PositionRadiansToDegrees(pos4)
#ra4 = ra4 * math.cos(math.radians(dec4))

PrintAndLog('RA1, Dec1: (%f, %f)' % (ra1, dec1), logfile)
PrintAndLog('RA2, Dec2: (%f, %f)' % (ra2, dec2), logfile)
PrintAndLog('RA3, Dec3: (%f, %f)' % (ra3, dec3), logfile)
PrintAndLog('RA4, Dec4: (%f, %f)' % (ra4, dec4), logfile)

minRA = min(ra1, ra2, ra3, ra4)
minDec = min(dec1, dec2, dec3, dec4)
maxRA = max(ra1, ra2, ra3, ra4)
maxDec = max(dec1, dec2, dec3, dec4)

#centerDec = ((maxDec - minDec) / 2.0) + minDec
#minRA = minRA * math.cos(math.radians(centerDec))
#maxRA = maxRA * math.cos(math.radians(centerDec))

# Because some CCDS are rotated by the WCS, make sure the min is actually the min.
if maxRA < minRA:
	temp = minRA
	minRA = maxRA
	maxRA = temp
if maxDec < minDec:
	temp = minDec
	minDec = maxDec
	maxDec = temp

PrintAndLog('Image limits', logfile)
PrintAndLog('(%f, %f):(%f, %f)' % (minRA, minDec, maxRA, maxDec), logfile)
if variables['useObservingList'] == 1:
	PrintAndLog('Using WCS from %s' % headerFile, logfile)

imageAreaDegrees = (maxRA - minRA) * (maxDec - minDec)

# Remove all of the comments from the CFIS header, since these were causing problems when writing the header out using GalSim.
del gsHeader['COMMENT']
del gsHeader['HISTORY']

gsHeader['seeing'] = imageSeeing

cfisImage.close()
            
numCopiesOfImage = int(variables['numCopiesOfImage'])
if numCopiesOfImage > 1:
	shears = GenerateShearValues(variables['e1s'], variables['e2s'])
else:
	shears = [[variables['e1s'], variables['e2s']]]

fullImages = []

for i in range(numCopiesOfImage):
	fullImages.append(galsim.ImageF(variables['image_size_x'], variables['image_size_y']))

bcoImage = None
#if int(variables['write_bco_file']) == 1:
#	bcoImage = galsim.ImageF(variables['image_size_x'], variables['image_size_y'])

psfImage = None
if int(variables['write_psf_file']) == 1:
	psfImage = galsim.ImageF(variables['image_size_x'], variables['image_size_y'])
	
missingPSFs = 0
galCount = 0
binCount = 0

truthFile = open('../Runs/%s/Image-Catalog-%i-%i-%i.txt' % (runName, coreNum, imageNum, ccdNum), 'w')
truthFile.write('Gal#, X, Y, r-Mag, HLR,  Sersic Index, e1 (intrinsic), e2 (intrinsic) \n')

# Generate the sheared PSF.
psfHLR = imageSeeing
psf = galsim.Moffat(beta = 3, fwhm = psfHLR)
numStars = 0

numLinesRead = 0
# Read in parameters from input catalog.
with open(catalogFilename) as catalog:
	for line in catalog:
		numLinesRead += 1

		# Skip the first line.
		if '#' in line:
			continue 

		# Create and add the galaxies to the image.
		params = line.split()
		dec = float(params[2])
		ra = float(params[1])# * math.cos(math.radians(dec))

		# Check to see if this object is within the bounds of the observation.
		if (ra >= minRA - OVERLAP_CHECK_DEG  and ra <= maxRA + OVERLAP_CHECK_DEG and dec >= minDec - OVERLAP_CHECK_DEG and dec <= maxDec + OVERLAP_CHECK_DEG):
			# It is, so get the rest of the parameters.
			mag = float(params[3])
			hlr = float(params[4])
			sInd = float(params[5])
			e1 = float(params[6])
			e2 = float(params[7])
		
			# Determine if the object is a galaxy or star.
			# TODO: Generalize this so that it can use flags from CFIS catalogs to determine whether galaxy or star.
			note = 'Galaxy'

			if sInd == -1:
				numStars += 1
				note = 'Star'

			AddGalaxyToImage(ra, minRA, dec, minDec, mag, hlr, sInd, e1, e2, galCount, psf, fullImages, psfImage, truthFile, note, bcoImage, wcs)
			galCount += 1

PrintAndLog('Number of lines in the complete catalog read: %i' % numLinesRead, logfile)
truthFile.close()

print('%i stars in this image.' % numStars)

noise = galsim.PoissonNoise(rng, sky_level=skyNoiseLevel)
noiseImage = galsim.ImageF(variables['image_size_x'], variables['image_size_y'])
noiseImage.addNoise(noise)

for i in range(numCopiesOfImage):
        # These images should all have identical noise, so add the noiseImage (which is a single realization of the noise model) to each of these.
        fullImages[i] += noiseImage
        fullImages[i] += int(variables['cfis_background'])
	#fullImages[i].write('../Runs/%s/FullImage-%i-%i-%iTEST.fits' % (runName, coreNum, imageNum, i))

	# Add the blank calibration rows and columns to the image.
	calColumns = galsim.ImageF(variables['calibration_columns'], variables['image_size_y'])
	calColumns.setOrigin(fullImages[i].origin)
	calColumns.setZero()
	bounds = calColumns.bounds & fullImages[i].bounds
        fullImages[i][bounds] = calColumns[bounds]


	calRows = galsim.ImageF(variables['image_size_x'], variables['calibration_columns'])
	calRows.setOrigin((1.0, variables['image_size_y'] - variables['calibration_columns'] + 1.0))
	calRows.setZero()
	bounds = calRows.bounds & fullImages[i].bounds
        fullImages[i][bounds] = calRows[bounds]

        fullImages[i].header = gsHeader
	fullImages[i].wcs = wcs
	fullImages[i].write('../Runs/%s/FullImage-%i-%i-%i-%i.fits' % (runName, coreNum, imageNum, i, ccdNum))

	if variables['write_noiseless_file'] == 1:
		noiseless = galsim.ImageF(variables['image_size_x'], variables['image_size_y'])
		noiseless = fullImages[i]
		noiseless -= noiseImage
        	noiseless.write('../Runs/%s/FullImageNoNoise-%i-%i-%i-%i.fits' % (runName, coreNum, imageNum, i, ccdNum))

        #fullImages[i] -= noiseImage
        #fullImages[i].write('../Runs/%s/FullImageNoNoise-%i-%i-%i.fits' % (runName, coreNum, imageNum, i))
        shearFile = open('../Runs/%s/FullImage-%i-%i-%i-%i-ShearValues.txt' % (runName, coreNum, imageNum, i, ccdNum), 'w')
        shearFile.write('e1 = %f, e2 = %f' % (shears[i][0], shears[i][1]))

if int(variables['write_psf_file']) == 1:
	psfImage.write('../Runs/%s/PSFImage-%i-%i-%i.fits' % (runName, coreNum, imageNum, ccdNum))

#if int(variables['write_bco_file']) == 1:
#	bcoImage.write('../Runs/%s/BCOImage.fits' % (runName))

if int(variables['write_weight_map']) == 1:
	weightMap = galsim.ImageF(variables['image_size_x'], variables['image_size_y'])
	weightMap.setZero()
	# TODO: Right now, weight maps are set to 1.  Need to calculate the inverse variance.  This might require normalised flat fields for the mosaic?
	weightMap += 1.0
	weightMap.write('../Runs/%s/FullImage-%i-%i-%i-%i-Weight.fits' % (runName, coreNum, imageNum, i, ccdNum))

if int(variables['write_flags_map']) == 1:
	flagsMap = galsim.ImageF(variables['image_size_x'], variables['image_size_y'])
	flagsMap.setZero()
	flagsMap.write('../Runs/%s/FullImage-%i-%i-%i-%i-Flags.fits' % (runName, coreNum, imageNum, i, ccdNum))

logfile.close()
