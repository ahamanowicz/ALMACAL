import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15 as cosmo
import astropy.units as u
from astropy import constants as const
def E(z):
	#using Planck15 cosmology and k=0 (flat universe)
	Om_M = cosmo.Om0 #Omega matter
	Om_L = 1-Om_M # omega Lambda , since Omega_k = 0
	#Om_k = 1- Om_M - Om_L # = 0 flat Universe

	E = np.sqrt(Om_M*(1+z)**3 + Om_L)
	return E

c =  const.c.to('km/s')
H0 =  cosmo.H(0)
DH = c/H0
print DH
z1,z2 = 0.1,0.11
d1,d2 = cosmo.luminosity_distance(z1),cosmo.luminosity_distance(z2)
print d1,d2, d2-d1

z1,z2 = 1.0 ,1.01
d1,d2 = cosmo.luminosity_distance(z1),cosmo.luminosity_distance(z2)
print d1,d2, d2-d1

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
rarcsec = 23.4796597742
deltaz=0.01
Z = np.arange(0.,3.,deltaz)
"""
fig1=plt.figure(1,figsize=(8,8))
plt.plot(Z,cosmo.luminosity_distance(Z), 'g', label='Luminocity distance')
plt.plot(Z, cosmo.comoving_distance(Z), 'b', label='Comving distance')
plt.plot(Z, cosmo.angular_diameter_distance(Z), 'r', label='Angular distance')
plt.xlabel("z")
plt.legend()
plt.ylabel("Mpc")
fig1.savefig("distances.pdf")
"""
#z1,z2=0.5,0.8
Vc = np.array([])
Vlum,Vang,Vc2, VN = np.array([]),np.array([]),np.array([]),np.array([])
deltaz=0.01

RL,RC,RA=np.array([]),np.array([]),np.array([])
print cosmo.kpc_proper_per_arcmin(rarcsec/60.)

for z in np.arange(0,3.,deltaz):
	z1 = z
	z2 = z1 + deltaz
	dL1,dL2 = cosmo.luminosity_distance(z1),cosmo.luminosity_distance(z2) #lum
	dC1,dC2 = cosmo.comoving_distance(z1), cosmo.comoving_distance(z2) # comoving
	dA1, dA2 = cosmo.angular_diameter_distance(z1), cosmo.angular_diameter_distance(z2) #angular

	#using angular distance for physical R
	#R1  = r * dA1
	#R1= r * dA2


	RL1 = (r * dL1/ (1.+z1)**2 ) #Mpc angular 
	RL2 = (r * dL2/ (1.+z2)**2 ) #Mpc angular
	RA1 = r*dA1
	RC1 = r *dC1
	
	RL = np.append(RL,RL1)
	RA = np.append(RA,RA1)
	RC = np.append(RC,RC1)

	#R_N
	z_mean = (z2+z1)/2.
	print z_mean
	p = cosmo.kpc_proper_per_arcmin(z_mean) #conversion factor from arcsmnin to kpc
	#angular size to physical
	R = p*(rarcsec/60.)*u.arcmin #kpc 
	# to comoving
	Rc = (R*(1+z_mean)).to(u.Mpc)
	Vol_lum = np.pi * abs(dL2-dL1) * RL2**2
	# angular volume
	Vol_ang = np.pi * abs(dA2-dA1) * RL2**2
	# comoving volume
	Vol_com = np.pi * abs(dC2-dC1) * RL2**2

	Vol_N = np.pi * abs(dC2-dC1)*Rc**2

	#Vc1 = np.pi * d1 * R2**2*(1+z1)**3 
	#Vc2 = np.pi * d2 * R2**2*(1+z2)**3 
	#Vc = Vc2-Vc1
	#print Vc1, Vc2
	#Vol_max = np.pi * abs(d2-d1) * R2**2
	
	Om_beam = np.pi * r**2
	
	element = 0
	nz = 1000.
	dz = (z2-z1)/nz
	print "o"
	
	for zz in np.linspace(z1,z2,nz):
		dL = cosmo.luminosity_distance(zz)

		element+= dL**2 / (1+zz)**2 / E(zz) * dz

	Vc_int = element*DH * Om_beam
	#Vc_other = DH * Om_beam* cosmo.luminosity_distance(z2)**2 / (1+z2)**2 / E(z2) * (z2-z1)
	
	Vc  = np.append(Vc, Vc_int)
	Vlum = np.append(Vlum, Vol_lum)
	Vang = np.append(Vang, Vol_ang)
	Vc2 = np.append(Vc2, Vol_com)
	VN = np.append(VN,Vol_N)
fig2 = plt.figure(2,figsize=(8,8))
Z = np.arange(0.,3.,deltaz)
plt.subplot(211)
plt.plot(Z, Vc, 'b', label = 'comoving integral')
plt.plot(Z, Vlum, 'g', label = 'luminocity')
#plt.plot(Z, Vang, 'r', label = 'angular diameter')
#plt.plot(Z, Vc2, 'm', label = 'comoving (distance)')
plt.plot(Z, VN, 'y', label = 'proper radius conversion')
plt.legend()
plt.xlabel("z")
plt.ylabel("Volume [Mpc^3]")
plt.subplot(212)
#plt.plot(Z, RL, 'g', label = 'luminocity')
#plt.plot(Z, RA, 'r', label = 'angular diameter', alpha=0.5)
#plt.plot(Z, RC, 'm', label = 'comoving ',alpha=0.5)
plt.plot(Z, Vc/Vlum, 'g', label= 'i = lum')
#plt.plot(Z, Vc/Vang, 'r', label= 'i = ang')
#plt.plot(Z, Vc/Vc2, 'm', label= 'i = com')
plt.plot(Z,Vc/VN, 'y',label='i = proper radius conv')
plt.xlabel("z")
plt.axhline(1,ls='--',c='gray')
plt.legend()
plt.ylabel("V_comoving / V_i")
plt.tight_layout()
fig2.savefig("volumes.pdf")
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