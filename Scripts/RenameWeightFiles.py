import glob
import os

def GetOutfileName(infile):
	folder = './' + infile.split('/')[1] + '/'
	return folder + infile.split('-')[-1]

def RenameFolderFiles(folder):
	weightFiles = glob.glob(folder + '*Weight.fits')

	for weightFile in weightFiles:
		os.rename(weightFile, weightFile.replace('sub-Weight', 'weight'))
		#print('%s, %s' % (weightFile, weightFile.replace('sub-Weight', 'weight')))

	weightFiles = glob.glob(folder + '*ShearValues.txt')

	for weightFile in weightFiles:
		os.rename(weightFile, weightFile.replace('sub-ShearValues', 'ShearValues'))
		#print('%s, %s' % (weightFile, weightFile.replace('sub-ShearValues', 'ShearValues')))

folders = ['0', '1', '2', '3']

for folder in folders:
	RenameFolderFiles('./' + folder + '/')
