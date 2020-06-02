from astropy.io import fits
import datetime
import galsim
import glob
import math
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy
import pyfits
import sys
from UtilityFunctions import ReadVariables

def PositionRadiansToDegrees(pos):
	ra = (pos.ra / galsim.radians) * (180.0 / numpy.pi)
	dec = (pos.dec / galsim.radians) * (180.0 / numpy.pi)
	return ra, dec

simFolder = sys.argv[1]

variables = ReadVariables()

simImages = glob.glob(simFolder + '*.fits')

minRA = 160.0
maxRA = 169.0
minDec = 29.0
maxDec = 31.5

for simImage in simImages:
	img = fits.open(simImage)
	header = img[0].header
	gsHeader = galsim.FitsHeader(header = header)
	wcs = galsim.wcs.readFromFitsHeader(gsHeader)[0]

	ra1, dec1 = PositionRadiansToDegrees(wcs.toWorld(galsim.PositionD(1.0, 1.0)))
	ra2, dec2 = PositionRadiansToDegrees(wcs.toWorld(galsim.PositionD(1.0, variables['image_size_y'])))
	ra3, dec3 = PositionRadiansToDegrees(wcs.toWorld(galsim.PositionD(variables['image_size_x'], variables['image_size_y'])))
	ra4, dec4 = PositionRadiansToDegrees(wcs.toWorld(galsim.PositionD(variables['image_size_x'], 1.0)))

	print('Rect((%f, %f), %f, %f' % (ra1, dec1, ra4 - ra1, dec1 - dec2))

	rect = patches.Rectangle((ra1, dec1), ra4 - ra1, dec1 - dec2, linewidth=1, edgecolor='r', facecolor='none')
	plt.gca().add_patch(rect)

plt.title('Chip locations')
plt.xlabel('RA')
plt.ylabel('Dec')
plt.xlim(minRA, maxRA)
plt.ylim(minDec, maxDec)
plt.gca().invert_xaxis()
plt.show()

