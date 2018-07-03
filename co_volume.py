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

def CO_Lum(z,fCO, Sline, dv=200):
	
	# z - line's redshift
	# fCO -  rest frequency of the chosen transition [GHz] 
	# Sline - strength of the signal n*rms [Jy]
	# dv - velocity span of the signal [km/s]

	# Lline [L_sol]

	dL = cosmo.luminosity_distance(z).value
	f_obs = fCO/(1+z)
	Lline = 1.04e-3 * Sline * dv * dL**2 * f_obs # L_sol
	return Lline*u.astrophys.solLum

def CO_Lum_prime(Lline, fCO):
	# LLine - result of CO_Lum 
	# fCO - rest frequency of the chosen transition [GHz] 

	Lline_prime = Lline/ 3.e-11 / fCO**3
	return Lline_prime
def CO_Lum_prime2(z,fCO, Sline, dv=200):
	#K km/s pc^2
	dL = cosmo.luminosity_distance(z).value
	f_obs = fCO/(1+z)
	Lline_prime = 3.25e7 * Sline * dv *  dL**2/(1+z)**3 / f_obs**2
	return Lline_prime

def sensitivity_function(r, mu, sigma):

	I = np.exp(-np.power(r - mu, 2.) / (2 * np.power(sigma, 2.)))

	return I
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

##### redshift coverage of the transition #######

#read in the Luminocity CO function values
Z,LCO = np.loadtxt('Luminocity_CO_prime_'+str(i)+'.txt', unpack=True,skiprows=2)
LCO *=Sline
co = CO[i] # CO transition central frequency [GHz]

z1,z2 = co/f2-1, co/f1 -1 #z1 lower, z2 higher redshift
dz = 4./100000.
# for this range of z calculate the range of LCO'
# safety
#LCO table
N = np.size(Z)

#full volume probled by the cube
nz = 100.
dd = (z2-z1)/nz
Vol_max=0
print r_max
dr = r_max/10.

d1 = cosmo.luminosity_distance(z1).value #dolna granica calkowania | fixed
d2 = cosmo.luminosity_distance(z2).value
print "f0",f0
for i in range(np.size(Z)):
	if abs(Z[i] - np.float(z2)) < 0.001 :
		Lmax= LCO[i]
	if abs(Z[i] - np.float(z1)) < 0.001 :
		Lmin = LCO[i]
print Lmax*1e-8, Lmin*1e-8
Vol_max = 0
print d1,d2
L_small = np.array([])
Z_small = np.array([])
for z in np.linspace(z1,z2,100):

	L_small=np.append(L_small,CO_Lum_prime2(z,co, Sline, dv=200))
	Z_small = np.append(Z_small,z)
#print L_small
for r in np.linspace(0,r_max-dr,100): # r in arcsec
	cr = sensitivity_function(r,mu=0.,sigma=sigma)#correction on intensity function rin arcsec
	L = Lmin * cr
	print L
	for i in range(np.size(Z_small)):
		#print L_small[i]
		if abs(L_small[i] - np.float(L)) < 1e8 :
			Lnew= L_small[i]
			zlim=Z_small[i]
	#print L/1e8, Lnew/1e8
	d2 = cosmo.luminosity_distance(zlim).value
	dd = (d2-d1)/nz
	#print d2
"""
	if d2 > d1:

		for d in np.linspace(d1,d2-dd,nz):
				
			#field of a ring in Mpc^2
			z = z_at_value(cosmo.luminosity_distance, d*u.Mpc)
			R1 = (np.radians(r/3600.) * d/ (1.+z)**2 ) #Mpc r " -> radians
			R2 = (np.radians((r+dr)/3600.) * d/ (1.+z)**2 ) #Mpc 

			Vol_max +=np.pi*(R2**2-R1**2)*dd
			#print R1,R2	
			#print sensitivity_function(r,mu=0.,sigma=sigma)
		
		#print d,d2

print "Vol5", Vol_max

for z in np.linspace(z1,z2-dd,nz):
	dS=0	
	d = cosmo.luminosity_distance(z).value
	for r in np.linspace(0,r_max-dr,100): # r in arcsec
		cr = 1#sensitivity_function(r,mu=0.,sigma=sigma)#correction on intensity function rin arcsec
		#field of a ring in Mpc^2
		R1 = (np.radians(r/3600.) * d/ (1.+z)**2 ) #Mpc r " -> radians
		R2 = (np.radians((r+dr)/3600.) * d/ (1.+z)**2 ) #Mpc 
		dS+=np.pi*(R2**2 - R1**2)*cr
		#print sensitivity_function(r,mu=0.,sigma=sigma)
	d1 = cosmo.luminosity_distance(z).value
	d2 = cosmo.luminosity_distance(z+dd).value
	Vol_max += dS * (d2-d1) #Mpc**3 physical
Vol_com_max = Vol_max *(1+z2)**3# / 1.e9 #Gpc

print z1,z2, Vol_max, Vol_com_max
d1,d2 = cosmo.luminosity_distance(z1).value,  cosmo.luminosity_distance(z2).value
#print fov(f0), imsize_arcsec 
r1,r2 =  np.radians(rmax/3600.) * d1/ (1.+z1)**2 , np.radians(rmax/3600.) * d2/ (1.+z2)**2 
print d1,r1,d2,r2

if z2 > 0 : #proceed
	if z1 < 0: z1 = 0 #if the transition lower limit dont fit the cube we start integrating at 0 
	for n in range(N):
		if Z[n] < z1:
			Vol = 0
		if Z[n] < z2 and Z[n] > z1:
			dd = (Z[n]-z1)/nz
			Vol = 0
			for z in np.linspace(z1,Z[n]-dd,nz):

				d = cosmo.luminosity_distance(z).value
				dS=0
				for r in np.linspace(0,r_max-dr,10): # r in arcsec
					cr = sensitivity_function(r,mu=0.,sigma=sigma)#correction on intensity function rin arcsec
					#field of a ring in Mpc^2
					R1 = (np.radians(r/3600.) * d/ (1.+z)**2 ) #Mpc r " -> radians
					R2 = (np.radians((r+dr)/3600.) * d/ (1.+z)**2 ) #Mpc 
					dS+=np.pi*(R2**2 - R1**2)/ cr
				d1 = cosmo.luminosity_distance(z).value
				d2 = cosmo.luminosity_distance(z+dd).value
				Vol += dS * (d2-d1) #Mpc**3 physical
			Vol_com = Vol *(1+Z[n])**3 # Mpc***3
			plt.plot(LCO[n], Vol_com,'ro')
		if Z[n] > z2:
			if LCO[n] < 1e8:
				#whole volume element
				Vol = Vol_com_max
				plt.plot(LCO[n], Vol_com_max,'bo')
else: 	
	print "Transition outside the cube frequency range!"
	Volume = 0	
plt.show()
"""

