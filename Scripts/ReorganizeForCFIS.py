# REORGANIZES THE OUTPUT FROM THE SIMULATION CODE INTO A FORMAT THAT IS PROCESSABLE BY LUDO'S PIPELINE.
# THIS SCRIPT IS OUTSIDE THE SCOPE OF WHAT I WANTED TO DO FOR THE SIMULATION CODE.
# IT THEREFORE EXISTS AS A SEPARATE SCRIPT, NOT INTENDED TO BE DISTRIBUTED WITH THE SIMULATION CODE.
# THIS SCRIPT IS INTENDED TO BE RUN ON LUDO'S MACHINE (CONSUS).

from astropy.io import fits
import glob
import numpy as np
import os

TILE_NAME = 'CFIS_147.9_40.5'
RUN_NUMBER = 1
RUN_NAME = 'run_CFIS_GALSIM_%i' % RUN_NUMBER

def RewriteHeaders():
	headers = glob.glob('./%s/r.MP9602/headers_V1.1.0A/*.head' % TILE_NAME)

	for header in headers:
		RewriteHeaderFile(header)

def RewriteHeaderFile(headerFilename):
	RemoveComments(headerFilename)
	RewritePVValues(headerFilename)

def RemoveComments(headerFilename):
	infile = open(headerFilename, 'r')
	lines = infile.readlines()
	infile.close()
	outfile = open(headerFilename, 'w')
	
	for line in lines:
		if not ('COMMENT' in line or 'HISTORY' in line):
			outfile.writeline(line)

	outfile.close()

def RewritePVValues(headerFilename):
	header = fits.Header.fromtextfile(headerFilename)
	header['PV1_0'] = 0.0
	header['PV1_1'] = 1.0
	header['PV1_2'] = 0.0
	header['PV1_3'] = 0.0
	header['PV1_4'] = 0.0
	header['PV1_5'] = 0.0
	header['PV1_6'] = 0.0
	header['PV1_7'] = 0.0
	header['PV1_8'] = 0.0
	header['PV1_9'] = 0.0
	header['PV1_10'] = 0.0
	header['PV2_0'] = 0.0
	header['PV2_1'] = 1.0
	header['PV2_2'] = 0.0
	header['PV2_3'] = 0.0
	header['PV2_4'] = 0.0
	header['PV2_5'] = 0.0
	header['PV2_6'] = 0.0
	header['PV2_7'] = 0.0
	header['PV2_8'] = 0.0
	header['PV2_9'] = 0.0
	header['PV2_10'] = 0.0
	header['RUN'] = RUN_NAME

	header.totextfile(headerFilename, endcard=True, overwrite=True)

def GetOutfileName(infile, folder):
	return folder + infile.split('-')[-1]

def RenameFolderFiles(folder):
	fitsFiles = glob.glob(folder + '*sub.fits')
	txtFiles = glob.glob(folder + 'FullImage*.txt')

	for fitsFile in fitsFiles:
		os.rename(fitsFile, GetOutfileName(fitsFile, folder))

	for txtFile in txtFiles:
		os.rename(txtFile, GetOutfileName(txtFile, folder))

def RenameImages():
	RenameFolderFiles('./%s/r.MP9602/single_V1.1.0A/' % TILE_NAME)

def RenameWeightFiles(folder):
	weightFiles = glob.glob(folder + '*Weight.fits')

	for weightFile in weightFiles:
		catFileBefore = '-'.join(weightFile.split('-')[:-3]).replace('FullImage', 'Image-Catalog') + '.txt'
		catFileAfter = '-'.join(weightFile.split('-')[0:1]).replace('FullImage', 'Image-Catalog-') + ''.join(weightFile.split('-')[-2].split('.')[:-1]) + '.txt'
		os.rename(catFileBefore, catFileAfter)
		os.rename(weightFile, weightFile.replace('sub-Weight', 'weight').split('-')[-1])

def RenameWeights():
	RenameWeightFiles('./%s/r.MP9602/single_V1.1.0A/' % TILE_NAME)

def RewriteFlagBitpix(folder):
	weightFilename = folder + '%s_r.MP9602.V1.1.0A.swarp.cut.flag.fits' % TILE_NAME 
	hdu = fits.open(weightFilename)
	hdu[0].data = np.zeros(hdu[0].data.shape, dtype=np.uint8)
	hdu[0].header['BITPIX'] = 8
	hdu.writeto(weightFilename, overwrite=True)
	hdu.close()

def RewriteWeightBitpix(folder):
	weightFilename = folder + '%s_r.MP9602.V1.1.0A.swarp.cut.weight.fits' % TILE_NAME 
	hdu = fits.open(weightFilename)
	hdu[0].data = np.ones(hdu[0].data.shape, dtype=np.uint8)
	hdu[0].header['BITPIX'] = 8
	hdu.writeto(weightFilename, overwrite=True)
	hdu.close()

def ZipFlagSumWeightImages():
	flagFile = glob.glob('./%s/r.MP9602/coadd_V1.1.0A/CFIS*.flag.fits' % TILE_NAME)[0]
	sumFile = glob.glob('./%s/r.MP9602/coadd_V1.1.0A/CFIS*.sum.fits' % TILE_NAME)[0]
	weightFile = glob.glob('./%s/r.MP9602/coadd_V1.1.0A/CFIS*.weight.fits' % TILE_NAME)[0]

	os.system('gzip %s' % flagFile)
	os.system('gzip %s' % sumFile)
	os.system('gzip %s' % weightFile)

def MakeDirectories():
	os.makedirs('/home/ispitzer/data/r.MP9602/')
	os.makedirs('/home/ispitzer/data/r.MP9602/%s/' % RUN_NAME)
	os.makedirs('/home/ispitzer/data/r.MP9602/%s/SCIENCE_r.MP9602/' % RUN_NAME)
	ps.makedirs('/home/ispitzer/data/r.MP9602/%s/WEIGHT/' % RUN_NAME)

MakeDirectories()
RewriteHeaders()
RenameImages()
RenameWeights()
RewriteFlagBitpix('./%s/r.MP9602/coadd_V1.1.0A/' % TILE_NAME)
RewriteWeightBitpix('./%s/r.MP9602/coadd_V1.1.0A/' % TILE_NAME)
ZipFlagSumWeightImages()
