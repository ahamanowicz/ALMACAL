import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import z_at_value
import astropy.units as u

def fov(freq):
	#frequency in in GHz
	freq*=1e9
	c = 3.e8
	d = 12.
	FOV = np.degrees(1.22*c/d/freq)*3600 #in arcsec
	return FOV

def CO_Lum(Sline, dv):
	# Sline = strength of the signal n*rms [Jy]
	# dv  [km/s]
	# Lline [L_sol]

	Lline = 1.04e-3 * Sline * dv * dL**2 * f
	return Lline

def CO_Lum_prime(Lline, fJ):
	# LLine - result of CO_Lum 
	# fJ - rest frequency of the line

	Lline_prime = Lline/ 3.e-11 / fJ**3
	return Lline_prime

######### Chose CO transition [GHz]
# 0-CO(1-0) 1-CO (2-1) 2-CO(3-2) 3-CO(4-3) 4-CO(5-4) 5-CO(6-5) 6-CO(7-6)

CO = np.array([115.27,230.538, 354.796, 461.041, 576.268, 691.473, 806.652])

# transition index 
i = 2 # CO(3-2)

####### read data from the cube *stats file
#single file (later change to whole sample or write short bash script)
fname = 'uid___A001_X1b1_X74.cube_1.J1924-2914_B7.stats'
with open(fname) as f:
    content = f.readlines()
content = [x.strip() for x in content] 

rms =  float(content[1].split()[1])
f1,f2 = float(content[7].split()[2]), float(content[7].split()[3]) #GHz
chan_width = float(content[8].split()[2]) #GHz
beam_radius = float(content[5].split()[3])*np.sqrt(2)*float(content[2].split()[4]) # imsize [pix] * sqrt(2) * cell [''/pix]

if (f1 > f2) : f1,f2 = f2,f1 # order asending
f0 = abs(f1+f2)/2. # central frequency of the cube

FOV  = fov(f0) # freq in GHz FOV of the detector ["]


##### redshift coverage of the transition #######

co = CO[2] # CO transition central frequency [GHz]

z1,z2 = co/f1-1, co/f2 -1

# safety
if z2 > 0 : #proceed
	if z1 < 0: z1 = 0 #if the transition lower limit dont fit the cube we start integrating at 0 

	#single cone 
	Vol = 0.
	dL1, dL2 = cosmo.luminosity_distance(z1).value, cosmo.luminosity_distance(z2).value # [Mpc], limits of the integration 
	nd = 100 # 100 elements for the integral
	dist = np.linspace(dL1,dL2,nz) 
	dd = abs(dL2-dL1)/nd #thinckes of the cone slice Mpc

	
	#print z_at_value(cosmo.luminosity_distance, dL1*u.Mpc)

	#perfect detector case
	for d in dist:  
	
		# volume element
		z = z_at_value(cosmo.luminosity_distance, d*u.Mpc)
		dS = np.pi*(FOV * d/ (1.+z)**2 )**2 
		Volume += dS * dd #Mpc**3 physical
		CO_Lum



"""

else: 
	print "Transition outside the cube frequency range!"
	Volume = 0


####### signal redshift limit  ######
n_sigma = 5 
dv = 200 # signal width in km/s
c = 3e5 # speed of ight in km/s


#### integrate over observation cone - physical FOV changes with redshift
Z = 2.0 # depth
dz = 0.001
distance = np.arange(0,Z,dz)
Volume = 0.

for i in range(np.size(distance)):
	z = distance[i]
	d =  cosmo.luminosity_distance(z) #Mpc

	dS = np.pi*(FOV * d/ (1.+z)**2 )**2 #FOV in radians, dL in Gpc
	d_dist = cosmo.luminosity_distance(z+dz) - cosmo.luminosity_distance(z)
	Volume += dS * d_dist #Gpc**3 physical
vol = Volume.value/1.e9	
vol_com = vol * (1+Z)**3 #Gpc**3 comoving
print vol,"Gpc^3",vol_com, "Gpc^3"
"""