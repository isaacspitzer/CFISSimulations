import glob

headFilenames = glob.glob('../../Data/Ludo/headers_V1.1.0A/*.head')

for headFilename in headFilenames:
	headFile = open(headFilename, 'r')
	filename = headFilename.split('/')[-1]
	newHeadFile = open('../Runs/2020-01-03-13:54/headers/' + filename, 'w')

	for line in headFile:
		newLine = line

		if line.startswith('PV'):
			paramName = line.split('=')[0]
			val = line.split('=')[1].split('/')[0]
			comment = line.split('/')[1]

			if ('PV1_1' in line or 'PV2_1' in line) and not ('PV1_10' in line or 'PV2_10' in line):
				newLine = paramName + ' = 1.0 / ' + comment
			else:
				newLine = paramName + ' = 0.0 / ' + comment

			newHeadFile.write(newLine)
		else:
			newHeadFile.write(line)

	headFile.close()
	newHeadFile.close()
			
