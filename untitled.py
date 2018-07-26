import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u
def E(z):
	Om_M = 0.3089
	Om_L = 0.6911
	Om_k = 1- Om_M - Om_L # = 0 flat Universe

	E = np.sqrt(Om_M*(1+z)**3 + Om_k*(1+z)**2 + Om_L)
	return E


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
	#print 2*np.pi*dR*r*h, Volume
	Volume +=2*np.pi*dR*r*h
	

#print  Vcylinder, Volume, V
### comoving volume check

z1,z2 = 0.629784442242, 0.63511300924

Vc1 = 28.073857221*(u.Mpc)**3
r =  np.radians(23.4796597742/3600.)

#z1,z2=0.5,0.8
Vc = np.array([])
Vphys = np.array([])
for z in np.arange(0,3.,0.01):
	z1 = z
	z2 = z1 + 0.01
	d1,d2 = cosmo.luminosity_distance(z1),cosmo.luminosity_distance(z2)
	R1 = (r * d1/ (1.+z1)**2 ) #Mpc
	R2 = (r * d2/ (1.+z2)**2 ) #Mpc
	#Vc1 = np.pi * d1 * R2**2*(1+z1)**3 
	#Vc2 = np.pi * d2 * R2**2*(1+z2)**3 
	#Vc = Vc2-Vc1
	#print Vc1, Vc2
	Vol_max = np.pi * abs(d2-d1) * R2**2
	Om_beam = np.pi * r**2
	DH = 3000./0.7 * u.Mpc
	element = 0
	nz = 100.
	dz = (z2-z1)/nz

	for zz in np.linspace(z1,z2,nz):
		dL = cosmo.luminosity_distance(zz)

		element+= dL**2 / (1+zz)**2 / E(zz) * dz

	Vc_int = element*DH * Om_beam
	#Vc_other = DH * Om_beam* cosmo.luminosity_distance(z2)**2 / (1+z2)**2 / E(z2) * (z2-z1)
	Vc  = np.append(Vc, Vc_int)
	Vphys = np.append(Vphys, Vol_max)

fig = plt.figure(1,figsize=(8,8))
Z = np.arange(0.,3.,0.01)
plt.subplot(211)
plt.plot(Z, Vc, 'b', label = 'comoving')
plt.plot(Z, Vphys, 'g', label = 'geom')
plt.legend()
plt.xlabel("z")
plt.ylabel("Volume [Mpc^3]")
plt.subplot(212)
plt.plot(Z, Vc/Vphys, 'k', label= 'alpha')
plt.xlabel("z")
plt.ylabel("V_comoving / V_geom")
plt.tight_layout()
fig.savefig("volumes.pdf")
plt.show()
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