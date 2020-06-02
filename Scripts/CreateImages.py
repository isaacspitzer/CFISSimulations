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
#import pyfits
import scipy.stats as stats
from UtilityFunctions import ReadVariables
import string

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

def RemoveCommentsFromHeaderFile(headerFilename):
	infile = open(headerFilename, 'r')
	lines = infile.readlines()
	infile.close()

	outfile = open(headerFilename, 'w')
	for line in lines:
		if not ('COMMENT' in line or 'HISTORY' in line):
			outline = line.split('/')[0] + '\n'
			outline = ''.join(i for i in outline if ord(i)<128)
			outfile.write(outline)

	outfile.close()


def CheckForPVValues(variables, gsHeader):
	setValue = False
	for i in range(2):
		for j in range(11):
			if 'PV%i_%i' % (i + 1, j) in variables:
				setValue = True
				break

	if setValue:
		for i in range(2):
			for j in range(11):
				if 'PV%i_%i' % (i + 1, j) in variables:
					gsHeader['PV%i_%i' % (i + 1, j)] = variables['PV%i_%i' % (i + 1, j)]
				else:
					gsHeader['PV%i_%i' % (i + 1, j)] = 0.0

def GetOutputDirectory(rotationNum, shearNum, runName, observationName, imagedir):
	outdir = imagedir.replace('<DATE_TIME>', '%s' % runName).replace('<SHEAR_NUM>', '%i' % shearNum).replace('<ROTATION_NUM>', '%i' % rotationNum)
	if observationName != 'None':
		outdir = outdir.replace('<OBS_NAME>', '%s' % observationName)

	return outdir

# Takes a set of model parameters, generates a galaxy and adds it to the provided image.
def AddGalaxyToImage(ra, minRA, dec, minDec, mag, hlr, sInd, e1, e2, galCount, psf, fullImages, psfImage, truthFiles, note, bcoImage, wcs, numRotations, rotationAngle):
	galaxy = None

	# Calculate the location of the galaxy in pixel space.
	worldPosition = galsim.CelestialCoord(ra=ra * galsim.degrees, dec=dec * galsim.degrees)
	
	try:
		imagePosition = wcs.toImage(worldPosition)
	except:
		#PrintAndLog("Error adding object to image: RA: %f, Dec %f" % (ra, dec), logfile)
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
	galaxyOriginal = galaxy.withFlux(counts)
	shape = numpy.matrix('%f %f' % (e1, e2))
        
	for j in range(numRotations):
		# Give the galaxy it's intrinsic shape.
		cos2Phi = numpy.cos(2.0 * rotationAngle * j)
		sin2Phi = numpy.sin(2.0 * rotationAngle * j)
		rotation = numpy.matrix('%f %f;%f %f' % (-1.0 * cos2Phi, -1.0 * sin2Phi, sin2Phi, -1.0 * cos2Phi))
		newShape = shape.dot(rotation)
		
		newE1 = newShape.item(0)
		newE2 = newShape.item(1)

		gal_shape = galsim.Shear(g1=newE1, g2=newE2)
		galaxy = galaxyOriginal.shear(gal_shape)

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
		
			bounds = stamp.bounds & fullImages[j][i].bounds

			try:
				fullImages[j][i][bounds] += stamp[bounds]
			except:
				# If the galaxy stamp is too far outside the frame, this will fail.  Just return None.
				#PrintAndLog('Galaxy too far out of frame.  Not adding.', logfile)
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
				if j == 0:
					truthFiles[j].write('%i %f %f %f %f %f %f %f %s %f %f\n' % (galCount, imagePosition.x, imagePosition.y, mag, hlr, sInd, e1, e2, note, ra, dec))

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

print('CreateImages Version')
print(sys.version)

# We don't want to just draw galaxies that are 'centered' in the frame.  We also want galaxies that are centered outside the frame, but which overlap with the borders of the frame and bleed into the image.  This setting dictates how far outisde of the frame (in degrees) a galaxy can be and still be added to the image.
OVERLAP_CHECK_DEG = 0.0015

# Specify an output folder for all of the generated images.
coreName = coreNum

variables = ReadVariables()

numRotations = int(variables['num_rotations'])
rotationAngle = variables['rotation_angle']
rotationAngle = rotationAngle * (numpy.pi / 180.0)

outfolder = variables['outdir'] #'../Runs/%s/Images/' % (runName)

# Create the folder for the images.
if not os.path.isdir(outfolder + '%s/' % runName):
	try:
		os.makedirs(outfolder + '%s/' % runName)
		print('Created %s/%s/' % (outfolder, runName))
	except:
		print('Output folder was created by another process.')

logfile = open(outfolder + '%s/log-CreateImages-%i-%i.txt' % (runName, imageNum, ccdNum), 'w')
msg = 'Creating output folders'
PrintAndLog(msg, logfile)

try:
	os.makedirs('../Runs/%s/Images/Catalogs/' % runName)
	os.makedirs('../Runs/%s/Headers/' % runName)
	
except:
	PrintAndLog('Catalog folder already created by another process.', logfile)

numCopiesOfImage = int(variables['numCopiesOfImage'])
if numCopiesOfImage > 1:
	shears = GenerateShearValues(variables['e1s'], variables['e2s'])
else:
	shears = [[variables['e1s'], variables['e2s']]]

imagedir = outfolder + variables['imagedir']

for j in range(numCopiesOfImage):
	for i in range(numRotations):
		observationName = variables['observationName']
		outdir = GetOutputDirectory(i, j, runName, observationName, imagedir)

		# Create the folder for the images.
		if not os.path.isdir(outdir):
			try:
				os.makedirs(outdir)
				os.makedirs(outdir.replace('single', 'headers'))
			except:
				print('Output folder was created by another process.')


# Make a local copy of the config file for future reference.
if not os.path.isfile(outfolder + '%s/Variables.cfg' % runName):
	copyfile('Variables.cfg', outfolder + '%s/Variables.cfg' % runName)

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
#elif variables['useObservingList'] == 2:
#	cfisImage = fits.open(variables['dataPath'] + headerFile.replace('.head', 'C.sub.fits'))
elif variables['useObservingList'] == 0:
	cfisImage = fits.open(variables['wcs_image'])

cfisHeader = None

print('Getting header')
if cfisImage != None:
	if int(variables['numCCDs']) > 1:
		cfisHeader = cfisImage[ccdNum+1].header
	else:
		cfisHeader = cfisImage[0].header

if 'headerPath' in variables:
	print('Reading from %s' % (variables['headerPath'] + headerFile.split('/')[-1].replace('C.sub.fits', '.head')))

	#headerTextFile = open(variables['headerPath'] + headerFile.split('/')[-1].replace('C.sub.fits', '.head'), 'r')

	if cfisHeader != None:
		cfisHeader.clear()

	#for line in headerTextFile:
	#	varName = line.split('=')[0].rstrip()
	#	if 'HISTORY' not in varName and 'COMMENT' not in varName and 'END' not in varName:
	#		try:
	#			cfisHeader[varName] = float(line.split('=')[1].split('/')[0])
	#		except:
	#			cfisHeader[varName] = line.split('=')[1].split('/')[0]
	#
	#headerTextFile.close()
	RemoveCommentsFromHeaderFile(variables['headerPath'] + headerFile.split('/')[-1].replace('C.sub.fits', '.head'))
	#newHead = fits.Header.fromtextfile(variables['headerPath'] + headerFile.split('/')[-1].replace('C.sub.fits', '_NoComments.head'), endcard=False, padding=False)
	newHead = fits.Header.fromtextfile(variables['headerPath'] + headerFile.split('/')[-1].replace('C.sub.fits', '.head'))
	cfisHeader = newHead

# Open the header and rewrite the telescope pointing.
#gsHeader = galsim.FitsHeader(header = cfisHeader)

# TODO: Approximating as flat rectangle.
#if ra >= 0.0:
#	c = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
#	cHMS = c.to_string('hmsdms')
#	raString = cHMS.split()[0].replace('h', ':').replace('m', ':').replace('s', '')
#	decString = cHMS.split()[1].replace('d', ':').replace('m', ':').replace('s', '').replace('+', '')
#	#gsHeader['ra'] = raString
#	cfisHeader['ra'] = raString
#	#gsHeader['dec'] = decString
#	cfisHeader['dec'] = decString
#	#gsHeader['crval1'] = ra
#	cfisHeader['crval1'] = ra
#	#gsHeader['crval2'] = dec
#	cfisHeader['crval2'] = dec
#	#gsHeader['ra_deg'] = ra
#	cfisHeader['ra_deg'] = ra
#	#gsHeader['dec_deg'] = dec
#	cfisHeader['dec_deg'] = dec

#gsHeader['datasec'] = '[%i:%i,1:%i]' % (variables['calibration_columns'] + 1, variables['image_size_x'] - variables['calibration_columns'], variables['image_size_y'] - variables['calibration_columns'])
cfisHeader['datasec'] = '[%i:%i,1:%i]' % (variables['calibration_columns'] + 1, variables['image_size_x'] - variables['calibration_columns'], variables['image_size_y'] - variables['calibration_columns'])
#gsHeader['fscale'] = variables['fscale']
cfisHeader['fscale'] = variables['fscale']

# Remove all of the comments from the CFIS header, since these were causing problems when writing the header out using GalSim.
try:
	#del gsHeader['COMMENT']
	del cfisHeader['COMMENT']
except:
	PrintAndLog('No COMMENTs to delete from header.', logfile)

try:
	#del gsHeader['HISTORY']
	del cfisHeader['HISTORY']
except:
	PrintAndLog('No HISTORYs to delete from header.', logfile)

cfisHeader['SEEING'] = imageSeeing
cfisHeader['GAIN'] = variables['gain']
cfisHeader['MAGZ'] = variables['zero_point_r']
#cfisHeader.totextfile(variables['headerPath'] + headerFile.split('/')[-1].replace('C.sub.fits', '.head'), endcard=True, overwrite=True)
cfisHeader.totextfile('../Runs/%s/Headers/' % (runName) + headerFile.split('/')[-1].replace('C.sub.fits', '.head'), endcard=True, overwrite=True)

if variables['overwrite_pv_values'] > 0.0:
	#CheckForPVValues(variables, gsHeader)
	CheckForPVValues(variables, cfisHeader)

gsHeader = galsim.FitsHeader(header = cfisHeader)
gsHeader['SEEING'] = imageSeeing
gsHeader['GAIN'] = variables['gain']
gsHeader['MAGZ'] = variables['zero_point_r']
#gsHeader['BSCALE'] = 1.0
#gsHeader['BZERO'] = 0.0
#gsHeader.write()

wcs = galsim.wcs.readFromFitsHeader(gsHeader)[0]
PrintAndLog("WCS is a %s" % wcs.__class__.__name__, logfile)
#PrintAndLog('gsHeader[PV1_1] == %s' % gsHeader['PV1_1'])

#  Want to find RA, Dec limits of the image.  Get the RA/Decs of each corner and find the min and max.
pos1 = wcs.toWorld(galsim.PositionD(1.0, 1.0))
pos2 = wcs.toWorld(galsim.PositionD(1.0, variables['image_size_y']))
pos3 = wcs.toWorld(galsim.PositionD(variables['image_size_x'], variables['image_size_y']))
pos4 = wcs.toWorld(galsim.PositionD(variables['image_size_x'], 1.0))
#print('POS1')
#print(pos1)
#print(pos2)
#print(pos3)
#print(pos4)
# Convert the radians to degrees.
ra1, dec1 = PositionRadiansToDegrees(pos1)
ra2, dec2 = PositionRadiansToDegrees(pos2)
ra3, dec3 = PositionRadiansToDegrees(pos3)
ra4, dec4 = PositionRadiansToDegrees(pos4)
#ra1 = pos1.x
#ra2 = pos2.x
#ra3 = pos3.x
#ra4 = pos4.x
#dec1 = pos1.y
#dec2 = pos2.y
#dec3 = pos3.y
#dec4 = pos4.y

#PrintAndLog('RA1, Dec1: (%f, %f)' % (ra1, dec1), logfile)
#PrintAndLog('RA2, Dec2: (%f, %f)' % (ra2, dec2), logfile)
#PrintAndLog('RA3, Dec3: (%f, %f)' % (ra3, dec3), logfile)
#PrintAndLog('RA4, Dec4: (%f, %f)' % (ra4, dec4), logfile)

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

if cfisImage != None:
	cfisImage.close()

fullImages = []

for j in range(numRotations):
	fullImageRotations = []
	for i in range(numCopiesOfImage):
		fullImageRotations.append(galsim.ImageF(variables['image_size_x'], variables['image_size_y']))

	fullImages.append(fullImageRotations)

bcoImage = None
#if int(variables['write_bco_file']) == 1:
#	bcoImage = galsim.ImageF(variables['image_size_x'], variables['image_size_y'])

psfImage = None
if int(variables['write_psf_file']) == 1:
	psfImage = galsim.ImageF(variables['image_size_x'], variables['image_size_y'])
	
missingPSFs = 0
galCount = 0
binCount = 0

truthFiles = []
inputFilename = headerFile.replace('./', '').replace('.fits', '')

#for i in range(numRotations):
if observationName != 'None':
	truthFile = open('../Runs/%s/Images/Catalogs/%s-Image-Catalog.txt' % (runName, inputFilename), 'w')
else:
	truthFile = open('../Runs/%s/Images/Catalogs/Image-Catalog-%i-%i-%i-%s.txt' % (runName, coreNum, imageNum, ccdNum, inputFilename), 'w')		
		
truthFile.write('Gal#, X, Y, r-Mag, HLR,  Sersic Index, e1 (intrinsic), e2 (intrinsic) \n')
truthFiles.append(truthFile)

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

			AddGalaxyToImage(ra, minRA, dec, minDec, mag, hlr, sInd, e1, e2, galCount, psf, fullImages, psfImage, truthFiles, note, bcoImage, wcs, numRotations, rotationAngle)
			galCount += 1

PrintAndLog('Number of lines in the complete catalog read: %i' % numLinesRead, logfile)
truthFile.close()

print('%i stars in this image.' % numStars)

noise = galsim.PoissonNoise(rng, sky_level=skyNoiseLevel)
noiseImage = galsim.ImageF(variables['image_size_x'], variables['image_size_y'])
noiseImage.addNoise(noise)

for j in range(numRotations):
	for i in range(numCopiesOfImage):
		# These images should all have identical noise, so add the noiseImage (which is a single realization of the noise model) to each of these.
		fullImages[j][i] += noiseImage
		fullImages[j][i] += int(variables['cfis_background'])
		#fullImages[j][i].write('../Runs/%s/FullImage-%i-%i-%iTEST.fits' % (runName, coreNum, imageNum, i))

		# Add the blank calibration rows and columns to the image.
		calColumns = galsim.ImageF(variables['calibration_columns'], variables['image_size_y'])
		calColumns.setOrigin(fullImages[j][i].origin)
		calColumns.setZero()
		bounds = calColumns.bounds & fullImages[j][i].bounds
		fullImages[j][i][bounds] = calColumns[bounds]


		calRows = galsim.ImageF(variables['image_size_x'], variables['calibration_columns'])
		calRows.setOrigin((1.0, variables['image_size_y'] - variables['calibration_columns'] + 1.0))
		calRows.setZero()
		bounds = calRows.bounds & fullImages[j][i].bounds
		fullImages[j][i][bounds] = calRows[bounds]

		#outfolder = '../Runs/%s/Images/Shear_%i/%s_Rotation_%i/' % (runName, j, observationName, i)
		fullImages[j][i].header = gsHeader
		fullImages[j][i].wcs = wcs

		outdir = GetOutputDirectory(j, i, runName, observationName, imagedir)

		fullImages[j][i].write(outdir + '%s.fits' % (inputFilename))
		#else:
		#	fullImages[j][i].write('../Runs/%s/Images/Shear_%i/Rotation_%i/FullImage-%i-%i-%i-%s.fits' % (runName, i, j, coreNum, imageNum, ccdNum, inputFilename))

		if variables['write_noiseless_file'] == 1:
			noiseless = galsim.ImageF(variables['image_size_x'], variables['image_size_y'])
			noiseless = fullImages[j][i]
			noiseless -= noiseImage
			#if observationName != 'None':
			noiseless.write(outdir + '%s.fits' % (inputFilename))
			#else:
			#	noiseless.write('../Runs/%s/Images/Shear_%i/Rotation_%i/FullImageNoNoise-%i-%i-%i-%s.fits' % (runName, i, j, coreNum, imageNum, ccdNum, inputFilename))

		#if observationName != 'None':
		shearFile = open(outdir + '%s-ShearValues.txt' % (inputFilename), 'w')
		#else:
		#	shearFile = open('../Runs/%s/Images/Shear_%i/Rotation_%i/FullImage-%i-%i-%i-%s-ShearValues.txt' % (runName, i, j, coreNum, imageNum, ccdNum, inputFilename), 'w')
		shearFile.write('e1 = %f, e2 = %f' % (shears[i][0], shears[i][1]))

		if int(variables['write_psf_file']) == 1:
			psfImage.write(outdir + 'PSFImage-%i-%i-%i-%s.fits' % (coreNum, imageNum, ccdNum, inputFilename))

		#if int(variables['write_bco_file']) == 1:
		#	bcoImage.write('../Runs/%s/BCOImage.fits' % (runName))

		if int(variables['write_weight_map']) == 1:
			weightMap = galsim.ImageF(variables['image_size_x'], variables['image_size_y'])
			weightMap.setZero()
			weightMap += 1.0

			#if observationName != 'None':
			weightMap.write(outdir + '%s-Weight.fits' % (inputFilename))
			#else:
			#	weightMap.write('../Runs/%s/Images/Shear_%i/Rotation_%i/FullImage-%i-%i-%i-%s-Weight.fits' % (runName, i, j, coreNum, imageNum, ccdNum, inputFilename))

		if int(variables['write_flags_map']) == 1:
			flagsMap = galsim.ImageF(variables['image_size_x'], variables['image_size_y'])
			flagsMap.setZero()

			#if observationName != 'None':
			flagsMap.write(outdir + '%s-Flags.fits' % (inputFilename))
			#else:
			#	flagsMap.write('../Runs/%s/Images/Shear_%i/Rotation_%i/FullImage-%i-%i-%i-%s-Flags.fits' % (runName, i, j, coreNum, imageNum, ccdNum, inputFilename))

		if int(variables['uniqueNoise']) > 0:		
			noise = galsim.PoissonNoise(rng, sky_level=skyNoiseLevel)
			noiseImage = galsim.ImageF(variables['image_size_x'], variables['image_size_y'])
			noiseImage.addNoise(noise)

logfile.close()
