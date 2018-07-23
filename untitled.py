import numpy as np
import matplotlib.pyplot as plt

# cone
d1=3.
d2=10.
#H=10.
R2 = 5.
R1=1.
V = 1./3. * np.pi*(d2-d1)*(R2**2 + R2*R1 + R1**2)
Vhole = 1./3.*np.pi * d2*R2**2
nz=100.
i=0.
dR=R2/nz
Volume=0

Vcylinder = np.pi * R2**2 *(d2-d1)

h=d2-d1
for r in np.linspace(0,R2,nz): # r in arcsec
	print 2*np.pi*dR*r*h, Volume
	Volume +=2*np.pi*dR*r*h
	

print  Vcylinder, Volume, V


"""
#z in np.linspace(z1,z2-dd,nz):
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
"""