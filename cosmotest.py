from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import z_at_value
import numpy as np
import astropy.units as u

def sensitivity_function(r, mu, sigma):

	I = np.exp(-np.power(r - mu, 2.) / (2 * np.power(sigma, 2.)))

	return I
def CO_Lum_prime2(z,fCO, Sline, dv=200):
	#K km/s pc^2
	dL = cosmo.luminosity_distance(z).value
	f_obs = fCO/(1+z)
	Lline_prime = 3.25e7 * Sline * dv *  dL**2/(1+z)**3 / f_obs**2
	return Lline_prime
r_max = 30.
sigma=15.
n=10
z1,z2=0.1,0.3
nz = 100
dd = (z2-z1)/nz
dr = r_max/nz

"""
d1 = cosmo.luminosity_distance(z1).value
d2 = cosmo.luminosity_distance(z2).value
R1 = (np.radians(r_max/3600.) * d1/ (1.+z1)**2 )
R2 = (np.radians(r_max/3600.) * d2/ (1.+z2)**2 ) 
Vol = 1./3.*np.pi*(d2-d1)*(R2**2 +R1*R2+R1**2)
print "Vol1", Vol



Vol_max=0
dS=0


#dd = (D2-D1)/nz

for z in np.linspace(z1,z2-dd,nz):
 	#z_at_value(cosmo.luminosity_distance,d*u.Mpc)
 	d=cosmo.luminosity_distance(z).value
 	
	R = (np.radians(r_max/3600.) * d/ (1.+z)**2 )

	d2 = cosmo.luminosity_distance(z+dd).value
	R2 = (np.radians(r_max/3600.) * d2/ (1.+z+dd)**2 )
	Vol_max+=np.pi*(R**2)*(d2-d)
	
print "Vol2", Vol_max
Vol_max_2 = 0
for z in np.linspace(z1,z2-dd,nz):
	dS=0	
	d = cosmo.luminosity_distance(z).value
	for r in np.linspace(0,r_max-dr,100): # r in arcsec
		cr = sensitivity_function(r,mu=0.,sigma=sigma)#correction on intensity function rin arcsec
		#field of a ring in Mpc^2

		R1 = (np.radians(r/3600.) * d/ (1.+z)**2 ) #Mpc r " -> radians
		R2 = (np.radians((r+dr)/3600.) * d/ (1.+z)**2 ) #Mpc 

		dS +=2*np.pi*R2*(R2-R1)*cr
		#print R1,R2	
		#print sensitivity_function(r,mu=0.,sigma=sigma)
	
	d2 = cosmo.luminosity_distance(z+dd).value
	#print d,d2
	Vol_max_2 += dS *(d2-d) #Mpc**3 physical
print "Vol3",Vol_max_2

Vol_max_3 = 0



d1 = cosmo.luminosity_distance(z1).value
d2 = cosmo.luminosity_distance(z2).value
dx=(d2+d1)/2
zx=(z2+z1)/2.
dR = (np.radians(dr/3600.) * dx/ (1.+zx)**2 )
for r in np.linspace(0,r_max-dr,nz): # r in arcsec
	cr = sensitivity_function(r,mu=0.,sigma=sigma)#correction on intensity function rin arcsec
	d1=cosmo.luminosity_distance(z1).value
	d2=cosmo.luminosity_distance(z2).value *cr
	h=d2-d1
	#z2 = z_at_value(cosmo.luminosity_distance,d2*u.Mpc)
	dR = (np.radians(dr/3600.) * d2/ (1.+z2)**2 )

	#field of a ring in Mpc^2

	R1 = (np.radians(r/3600.) * d1/ (1.+z1)**2 ) #Mpc r " -> radians
	R2 = (np.radians(r/3600.) * d2/ (1.+z2)**2 ) #Mpc 

	Vol_max_3 +=np.pi*np.sqrt(h**2 + (R2-R1)**2)*(R2+R1)*dR
	#print R1,R2	
	#print sensitivity_function(r,mu=0.,sigma=sigma)

	#print d2, d2*cr
	#print d,d2
print"Vol4", Vol_max_3
"""
CO = np.array([115.27,230.538, 354.796, 461.041, 576.268, 691.473, 806.652])

# transition index 
i = 2 # CO(3-2)
Z,LCO = np.loadtxt('Luminocity_CO_prime_'+str(i)+'.txt', unpack=True,skiprows=2)
Vol_max_4 = 0
nz = 100
d1 = cosmo.luminosity_distance(z1).value #dolna granica calkowania | fixed
d2 = cosmo.luminosity_distance(z2).value
f0 = 282.4237
for i in range(np.size(Z)):
	if abs(Z[i] - np.float(z2)) < 0.0001 :
		Lmax= LCO[i]
	if abs(Z[i] - np.float(z1)) < 0.0001 :
		Lmin = LCO[i]
print Lmax*1e-8, Lmin*1e-8

for r in np.linspace(0,r_max-dr,100): # r in arcsec
	cr = sensitivity_function(r,mu=0.,sigma=sigma)#correction on intensity function rin arcsec
	L = Lmin * cr
	for i in range(np.size(Z)):
		if abs(LCO[i] - np.float(L)) < 1e8 :
			Lnew= LCO[i]
			zlim=Z[i]
	d2 = cosmo.luminosity_distance(zlim).value
	dd = (d2-d1)/nz
	print d2

	if d2 > d1:

		for d in np.linspace(d1,d2-dd,nz):
				
			#field of a ring in Mpc^2
			z = z_at_value(cosmo.luminosity_distance, d*u.Mpc)
			R1 = (np.radians(r/3600.) * d/ (1.+z)**2 ) #Mpc r " -> radians
			R2 = (np.radians((r+dr)/3600.) * d/ (1.+z)**2 ) #Mpc 

			Vol_max_4 +=np.pi*(R2**2-R1**2)*dd
			#print R1,R2	
			#print sensitivity_function(r,mu=0.,sigma=sigma)
		
		#print d,d2

print "Vol5", Vol_max_4
