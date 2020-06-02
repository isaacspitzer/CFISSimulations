# Read in the variables defined in './Variables.cfg'.
# The variables will be stored in a dictionary and returned.
# The function must ignore blank lines and handle Python style commenting with the # character.
def ReadVariables():
	variables = dict()

	with open("./Variables.cfg") as f:
		for line in f:
			if line.startswith('#') or not line.strip():
				continue

			noncommentPart = line.split('#')[0]
			eqInd = noncommentPart.find('=')
			varName = noncommentPart[:eqInd].strip()
		
			try:
				varVal = float(noncommentPart[eqInd + 1:].strip())
			except:
				varVal = noncommentPart[eqInd + 1:].strip()

			variables[varName] = varVal

	return variables

