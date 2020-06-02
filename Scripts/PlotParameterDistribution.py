import matplotlib.pyplot as plt
import numpy as np
from UtilityFunctions import ReadVariables

# Read in the catalog of galaxy models.
# This function will return a list of the magnitudes of the models, as well as a list of Galaxy objetcs containing the complete set of model properties.
def ReadModelCatalog(filename):
	mags = []
	hlrs = []
	sinds = []
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
			hlrs.append(hlr)
			sinds.append(sInd)
                        #models.append(Galaxy(hlr, sInd, 0.0, 0, 0, mag, axisRatio))			

                line = gemsfile.readline()

        return (mags, hlrs, sinds)

def ReadSimCatalog(filename):
	mags = []
	hlrs = []
	sinds = []

	lineCount = 0
	catfile = open(filename, 'r')

	for line in catfile:
		if lineCount == 0:
			lineCount = lineCount + 1
			continue

		data = line.split()
		mag = float(data[3])
		hlr = float(data[4])
		sind = float(data[5])
		
		mags.append(mag)
		hlrs.append(hlr)
		sinds.append(sind)

	return (mags, hlrs, sinds)

variables = ReadVariables()
cosmosMags, cosmosHlrs, cosmosSinds = ReadModelCatalog('../COSMOS/asu.csv')
simMags, simHlrs, simSinds = ReadSimCatalog('../Runs/2019-07-08-14:51/SimulatedCatalog.cat')

# Convert to numpy arrays.
#cosmosMags = np.asarray(cosmosMags)
#cosmosHlrs = np.asarray(cosmosHlrs)
#cosmosSinds = np.asarray(cosmosSinds)
#simMags = np.asarray(simMags)
#simHlrs = np.asarray(simHlrs)
#simSinds = np.asarray(simSinds)

# Normalize the distributions.
#cosmosMags = cosmosMags / float(len(cosmosMags))
#cosmosHlrs = cosmosHlrs / float(len(cosmosHlrs))
#cosmosSinds = cosmosSinds / float(len(cosmosSinds))
#simMags = simMags / float(len(simMags))
#simHlrs = simHlrs / float(len(simHlrs))
#simSinds = simSinds / float(len(simSinds))

plt.hist(cosmosMags, 50, normed=1, facecolor='green', alpha=0.5, label='COSMOS')
plt.hist(simMags, 50, normed=1, facecolor='blue', alpha=0.5, label='Simulation')
plt.xlabel('r-band magnitude')
plt.legend(loc=0)
plt.gca().set_yscale("log")
plt.show()
#plt.savefig('magnitude_distribution.eps')
plt.clf()

plt.hist(cosmosHlrs, 50, normed=1, facecolor='green', alpha=0.5, label='COSMOS')
plt.hist(simHlrs, 50, normed=1, facecolor='blue', alpha=0.5, label='Simulation')
plt.xlabel('Half-light radius')
plt.legend(loc=1)
plt.gca().set_yscale("log")
plt.show()
#plt.savefig('hlr_distribution.eps')
plt.clf()

plt.hist(cosmosSinds, 50, normed=1, facecolor='green', alpha=0.5, label='COSMOS')
plt.hist(simSinds, 50, normed=1, facecolor='blue', alpha=0.5, label='Simulation')
plt.xlabel('Sersic index')
plt.legend(loc=1)
plt.gca().set_yscale("log")
plt.show()
#plt.savefig('sersic_distribution.eps')
plt.clf()