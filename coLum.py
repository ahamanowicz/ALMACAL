import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import z_at_value
def CO_Lum_prime2(z,fCO, Sline, dv=200):
	#K km/s pc^2
	dL = cosmo.luminosity_distance(z).value
	f_obs = fCO/(1+z)
	Lline_prime = 3.25e7 * Sline * dv *  dL**2/(1+z)**3 / f_obs**2
	return Lline_prime
fname = 'uid___A001_X1b1_X74.cube_1.J1924-2914_B7.stats'
with open(fname) as f:
    content = f.readlines()
content = [x.strip() for x in content] 

rms =  float(content[1].split()[1])
f1,f2 = float(content[7].split()[2]), float(content[7].split()[3]) #GHz
chan_width = float(content[8].split()[2]) #GHz
beam_radius = float(content[5].split()[3])*np.sqrt(2)*float(content[2].split()[4]) # imsize [pix] * sqrt(2) * cell [''/pix]
imsize_pix = float(content[5].split()[3])
cell = float(content[2].split()[4])

imsize_arcsec = imsize_pix*cell # arces
imsize_radians = np.radians(imsize_arcsec/3600.) #to go into the cosmology equation this must be in radians
if (f1 > f2) : f1,f2 = f2,f1 # order asending
f0 = abs(f1+f2)/2. # central frequency of the cube

sigma = imsize_arcsec/2.35/1.5 # arcsec
r_max = 1.7*sigma #arcsec
print imsize_radians,sigma
Sline = 5 *rms # line strength
dv = 200 # [km/s] width of the line	
CO = np.array([115.27,230.538, 354.796, 461.041, 576.268, 691.473, 806.652])

Z=np.linspace(0,10,10000)
plt.figure(figsize=(10,6))
plt.plot(Z,np.log10(CO_Lum_prime2(Z,CO[0], Sline, dv=200)), label="CO(1-0)")
plt.plot(Z,np.log10(CO_Lum_prime2(Z,CO[1], Sline, dv=200)), label="CO(2-1)")
plt.plot(Z,np.log10(CO_Lum_prime2(Z,CO[2], Sline, dv=200)), label="CO(3-2)")
plt.plot(Z,np.log10(CO_Lum_prime2(Z,CO[3], Sline, dv=200)), label="CO(4-3)")
plt.plot(Z,np.log10(CO_Lum_prime2(Z,CO[4], Sline, dv=200)), label="CO(5-4)")
plt.plot(Z,np.log10(CO_Lum_prime2(Z,CO[5], Sline, dv=200)), label="CO(6-5)")
plt.plot(Z,np.log10(CO_Lum_prime2(Z,CO[6], Sline, dv=200)), label="CO(7-6)")
plt.axvline(1.5,ls='--',c='r')
#print zlim
plt.ylabel('log10(LCO prime)')
plt.xlabel('z lim')
plt.legend()
plt.grid()
plt.show()