import glob
import os

def GetOutfileName(infile):
	folder = './' + infile.split('/')[1] + '/'
	return folder + infile.split('-')[-1]

def RenameFolderFiles(folder):
	fitsFiles = glob.glob(folder + '*.fits')
	txtFiles = glob.glob(folder + 'FullImage*.txt')

	for fitsFile in fitsFiles:
		os.rename(fitsFile, GetOutfileName(fitsFile))
		#print('%s, %s' % (fitsFile, GetOutfileName(fitsFile)))

	for txtFile in txtFiles:
		os.rename(txtFile, GetOutfileName(txtFile))
		#print('%s, %s' % (txtFile, GetOutfileName(txtFile)))

folders = ['0', '1', '2', '3']

for folder in folders:
	RenameFolderFiles('./' + folder + '/')
