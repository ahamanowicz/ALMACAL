from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import z_at_value
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
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


Vol_max=0
nz=100
CO = np.array([115.27,230.538, 354.796, 461.041, 576.268, 691.473, 806.652])

# transition index 
i = 2 # CO(3-2)
Z,LCO = np.loadtxt('Luminocity_CO_prime_'+str(i)+'.txt', unpack=True,skiprows=2)

dr = r_max/nz
print dd
dS = 0

Vol = 0

d1 = cosmo.luminosity_distance(z1).value
d2 = cosmo.luminosity_distance(z2).value
print d1
dL_lim_max =0
z2=1.4
a = cosmo.luminosity_distance(1).value**2/(1+1)

print cosmo.luminosity_distance(1.5).value**2/(1+1.5), cosmo.luminosity_distance(2).value**2/(1+2)
#parabola
X1 = np.linspace(0,1.5,500)
Y1 = cosmo.luminosity_distance(X1).value**2/(1+X1)
coef1 = np.polyfit(X1,Y1,6)
#print coef1
p1 = np.poly1d(coef1)
#linear
X2 = np.linspace(1.5,8,1000)
Y2 = cosmo.luminosity_distance(X2).value**2/(1+X2)
coef2 = np.polyfit(X2,Y2,6)
print coef2
p2 = np.poly1d(coef2)
Z = np.linspace(0,8,1000)

y = 5.e7
print (p1-y).roots

for z in Z :
	func = cosmo.luminosity_distance(z).value**2/(1+z)
	plt.plot(z,func, 'ob', markersize=2, alpha=0.5)

plt.plot(Z,p1(Z), 'r', label='approx1')	
plt.plot(Z,p2(Z), 'g', label='approx2')	
plt.legend()
plt.xlabel("z")
plt.ylabel(r"$d_{L}^{2}(1+z)^{-1}$")
plt.show()

for i in range(np.size(Z)):
		if abs(Z[i] - np.float(z2)) < 1e-4:
			print Z[i], LCO[i] 

i=0
ddd = d2/nz
for r in np.linspace(0,r_max,nz): # r in arcsec
	cr = sensitivity_function(r, mu=0, sigma=sigma)
	dlim = d2*cr #upper limit
	d1_prime = ddd * i #lower limit
	if d1_prime < d1:
		d1_prime = d1
	
	h  = dlim - d1_prime # height of the element cilinder
	print dlim
	zlim = z_at_value(cosmo.luminosity_distance, dlim*u.Mpc)
	if h > 0:
		
		R = (np.radians(r/3600.) * dlim/ (1.+zlim)**2 )
		dR = (np.radians(dr/3600.) * dlim/ (1.+zlim)**2 )
		
		Vol += 2*np.pi*dR*R * h


	i+=1
print Vol
"""
for r in np.linspace(0,r_max,nz): # r in arcsec
	cr = sensitivity_function(r, mu=0, sigma=sigma)
	dlim = d2*cr
	dd = dlim - d1
	

	zlim = z_at_value(cosmo.luminosity_distance, dlim*u.Mpc)
	if dd > 0:
		print dlim
		R = (np.radians(r/3600.) * dlim/ (1.+zlim)**2 )
		dR = (np.radians(dr/3600.) * dlim/ (1.+zlim)**2 )
		
		Vol += 2*np.pi*dR*R * dd

Rmax = 	(np.radians(r_max/3600.) * d2/ (1.+z2)**2 )

Vol_real = Rmax**2 * d2 *np.pi/3.
print "shells",Vol, "Real", Vol_real

#dd = (D2-D1)/nz

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
"""