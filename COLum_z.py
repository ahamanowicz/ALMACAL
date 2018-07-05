#CO Luminocity function table calculaton

import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import z_at_value
def CO_Lum_prime(z,fCO, Sline, dv=200):
	#K km/s pc^2
	dL = cosmo.luminosity_distance(z).value
	f_obs = fCO/(1+z)
	Lline_prime = 3.25e7 * Sline * dv *  dL**2/(1+z)**3 / f_obs**2
	return Lline_prime


CO = np.array([115.27,230.538, 354.796, 461.041, 576.268, 691.473, 806.652])

 # CO transition central frequency [GHz
z1,z2=0,4
z_space = np.linspace(z1,z2,10000)
Sline = 1

for i in range(7):
	co = CO[i]
	g = open('Luminocity_CO_prime_'+str(i)+'.txt','w')
	g.write('#redshift \t LCOprime \n')
	for z in z_space:

		lum = CO_Lum_prime(z,co, Sline, dv=200)/1e8
	#	plt.plot(z,lum, 'ko')
		g.write(str(round(z,5))+'\t'+str(lum)+'\n')

	#plt.show()
	g.close()