import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import z_at_value
import astropy.units as u
########################################################
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

def L_CO_prime_inv(LCO,fCO,cr,Sline,dv=200):
	func = LCO*cr*fCO**2/(3.25e7 * Sline * dv )
 	return func

def sensitivity_function(r, mu, sigma):

	I = np.exp(-np.power(r - mu, 2.) / (2 * np.power(sigma, 2.)))

	return I
########################################################
### distance function:

# z <  1.5
X1 = np.linspace(0,1.5,1000)
Y1 = cosmo.luminosity_distance(X1).value**2/(1+X1)
coef1 = np.polyfit(X1,Y1,6)
p1 = np.poly1d(coef1)

# z > 1.5
X2 = np.linspace(1.5,5,2000)
Y2 = cosmo.luminosity_distance(X2).value**2/(1+X2)
coef2 = np.polyfit(X2,Y2,6)

p2 = np.poly1d(coef2)
Z = np.linspace(0,5,500)
N = np.size(Z)

##plot the distance function approximations
"""
plt.figure(1)
plt.plot(Z,np.log10(p1(Z)), 'y')
plt.plot(Z,np.log10(p2(Z)),'g')
plt.plot(Z,np.log10(cosmo.luminosity_distance(Z).value**2/(1+Z)),'r',alpha=0.5)
plt.xlabel("z")
plt.ylabel("log10(f)")
plt.show()
"""
######### Chose CO transition [GHz]
# 0-CO(1-0) 1-CO (2-1) 2-CO(3-2) 3-CO(4-3) 4-CO(5-4) 5-CO(6-5) 6-CO(7-6)

CO = np.array([115.27,230.538, 354.796, 461.041, 576.268, 691.473, 806.652])

# transition index 
i = 3 # CO(3-2)

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
#print imsize_radians,sigma
Sline = 5 *rms # line strength
dv = 200 # [km/s] width of the line 

##### redshift coverage of the transition  [cube limits] #######

co = CO[i] # CO transition central frequency [GHz]
z1,z2 = co/f2-1, co/f1 -1 #z1 lower, z2 higher redshift | limiting redshft from the band coverage
#z1,z2=1.7,2.2
logLUM = np.linspace(5,13,50)#np.array([1e6,1e7,1e8,1e9,1e10,1e11,1e12])
LUM = 10**logLUM
# security 
if z2 > 0:

	if z1 < 0 : z1=0
	d1 = cosmo.luminosity_distance(z1) # lower distance limit [Mpc]
	d2 = cosmo.luminosity_distance(z2) #higher distance limit [Mpc]

	#### INTEGRATION ####

	nz = 100.
	dd = (z2-z1)/nz
	print "r_max",r_max
	dr = r_max/nz

	## maximum value of the volume - turnicated cone
	R1 = (np.radians(r_max/3600.) * d1/ (1.+z1)**2 ) #Mpc
	R2 = (np.radians(r_max/3600.) * d2/ (1.+z2)**2 ) #Mpc
	Vol_max = 1./3. *np.pi * abs(d2-d1) *(R2**2 + R2*R1 + R1**2)
	Vol_max_2 = np.pi * abs(d2-d1) * R2**2
	print "R1,R2,d1,d2",R1,R2,d1,d2
	print "stoz",Vol_max, "walec", Vol_max_2
	
	#chosin the regime

	LCO_lim = CO_Lum_prime2(1.5,co, Sline, dv=200)
	print "LCOlim", np.log10(LCO_lim), z1,z2
	print d1,d2

	for l in LUM:
		Volume = 0.
		print Volume
		if  l < LCO_lim: 
			p = p1
			print "p1"
			###limiting redshifts from sensitivity
			#max
			cr=1
			func_z = L_CO_prime_inv(l,co,cr,Sline,dv=200)
			roots  = np.real((p-func_z).roots) #find solutions  and choose right root
			positive_roots =[ n for n in roots if n > 0]
			zlim_max = min(positive_roots) #in the center of an image

			#min
			cr = sensitivity_function(r_max, mu=0, sigma=sigma)
			func_z = L_CO_prime_inv(l,co,cr,Sline,dv=200)
			roots  = np.real((p-func_z).roots) #find solutions  and choose right root
			positive_roots =[ n for n in roots if n > 0]
			zlim_min = min(positive_roots) #at the edge of he image
		
			if (zlim_min > zlim_max) : zlim_min,zlim_max = zlim_max,zlim_min
			dlim_max = cosmo.luminosity_distance(zlim_max)
			dlim_min = cosmo.luminosity_distance(zlim_min)
			#print "zmax, zmin, LCO", zlim_max,zlim_min, np.log10(l)

		else: 
			p = p2
			print "p2"
			cr=1
			func_z = L_CO_prime_inv(l,co,cr,Sline,dv=200)
			roots  = np.real((p-func_z).roots) #find solutions  and choose right root
			positive_roots =[ n for n in roots if n > 1.5]
			zlim_max = min(positive_roots)

			#min
			cr = sensitivity_function(r_max, mu=0, sigma=sigma)
			func_z = L_CO_prime_inv(l,co,cr,Sline,dv=200)
			roots  = np.real((p-func_z).roots) #find solutions  and choose right root
			positive_roots =[ n for n in roots if n > 1.5]
			zlim_min = min(positive_roots)
		
			if (zlim_min > zlim_max) : zlim_min,zlim_max = zlim_max,zlim_min
			dlim_max = cosmo.luminosity_distance(zlim_max)
			dlim_min = cosmo.luminosity_distance(zlim_min)
		print "zmax, zmin, LCO", zlim_max,zlim_min, np.log10(l)
		print "dmax, dmin", dlim_max,dlim_min
		### CASE 1 - zlim_max > z2 and zlim_min  > z2; volume -> volume of the truncated cone
		if zlim_max > z2 and zlim_min  > z2:
			print "Case 1"
			Volume = Vol_max_2

		### CASE 2 - zlim_max < z1
		elif zlim_max < z1 :
			print "Case 2"
			Volume = 0
		### CASE 3 - zlim_max > z2 and zlim_min < z2  (always check if zlim_min > z1)

		elif zlim_max > z2 and zlim_min < z2:	
			print "Case 3"
			for r in np.linspace(0,r_max,nz): # r in arcsec
				cr = sensitivity_function(r, mu=0, sigma=sigma)
				func_z = L_CO_prime_inv(l,co,cr,Sline,dv=200)
		
				if  l < LCO_lim: 
					p = p1
					roots  = np.real((p-func_z).roots) #find solutions  and choose right root
					positive_roots =[ n for n in roots if n > 0]
				else  : 
					p = p2
					roots  = np.real((p-func_z).roots) #find solutions  and choose right root
					positive_roots =[ n for n in roots if n > 1.5]
				#print np.log10(l*cr),min(positive_roots)

				zlim=min(positive_roots)
				dlim = cosmo.luminosity_distance(zlim) #lower limit, defined by sensitivity
				
				if zlim  > z2:

					h=d2-d1
					R = (np.radians(r/3600.) * d2/ (1.+z2)**2 )
					dR = (np.radians(dr/3600.) * d2/ (1.+z2)**2 )
					Volume += 2*np.pi*dR*R*h
					
					
				elif zlim > z1:
				

					h = dlim - d1
					R = (np.radians(r/3600.) * dlim/ (1.+zlim)**2 )
					dR = (np.radians(dr/3600.) * dlim/ (1.+zlim)**2 )
					Volume += 2*np.pi*dR*R*h
				else:
					Volume +=0
				
					

		### CASE 4 - zlim_max < z2  and zlim_min > z1  (always check if zlim_min > z1)
		elif zlim_max < z2 and zlim_min  > z1:
			print "Case 4"
			#i=0
			## integrating by hollow cilinders

			for r in np.linspace(0,r_max,nz): # r in arcsec
				cr = sensitivity_function(r, mu=0, sigma=sigma)
				func_z = L_CO_prime_inv(l,co,cr,Sline,dv=200)
		
				if  l < LCO_lim: 
					p = p1
					roots  = np.real((p-func_z).roots) #find solutions  and choose right root
					positive_roots =[ n for n in roots if n > 0]
				else  : 
					p = p2
					roots  = np.real((p-func_z).roots) #find solutions  and choose right root
					positive_roots =[ n for n in roots if n > 1.5]
				#print np.log10(l*cr),min(positive_roots)

				zlim=min(positive_roots)
				dlim = cosmo.luminosity_distance(zlim) #lower limit, defined by sensitivity
				#print zlim
				h = dlim - d1
				R = (np.radians(r/3600.) * dlim/ (1.+zlim)**2 )
				dR = (np.radians(dr/3600.) * dlim/ (1.+zlim)**2 )
				Volume += 2*np.pi*dR*R*h
				
				#Volume -=1./3.*np.pi*(np.radians(r_max/3600.) * dlim_min/ (1.+zlim_min))**2*dlim_min 
		### CASE 5 - zlim_max < z2  and zlim_min < z1  (always check if zlim_min > z1)
		elif  zlim_max < z2 and zlim_max > z1 and zlim_min  < z1:
			print "Case 5"
			## integrating by hollow cilinders

			for r in np.linspace(0,r_max,nz): # r in arcsec
				cr = sensitivity_function(r, mu=0, sigma=sigma)
				func_z = L_CO_prime_inv(l,co,cr,Sline,dv=200)
				if  l < LCO_lim: 
					p = p1
					roots  = np.real((p-func_z).roots) #find solutions  and choose right root
					positive_roots =[ n for n in roots if n > 0]
				else  : 
					p = p2
					roots  = np.real((p-func_z).roots) #find solutions  and choose right root
					positive_roots =[ n for n in roots if n > 1.5]
				
				#print np.log10(l*cr),min(positive_roots)
				zlim=min(positive_roots)
				dlim = cosmo.luminosity_distance(zlim)
				if zlim > z1:
					
					#d1_prime = dlim_max/nz * i #lower limit
					h  = dlim - d1 # height of the element cilinder
					R = (np.radians(r/3600.) * dlim/ (1.+zlim)**2 )
					dR = (np.radians(dr/3600.) * dlim/ (1.+zlim)**2 )
					Volume += 2*np.pi*dR*R*h
					#print d1_prime, dlim
				else: 
					Volume += 0
				
		if Volume > Vol_max_2: print "WRONG"
		
		print "l,vol",np.log10(l), Volume
	
		plt.plot(np.log10(l), Volume, 'ob')
else: Volume = 0

plt.xlabel("log10(LCO_prime)")
plt.ylabel("Volume [Mpc^3]")
plt.show()


"""
else:
			if zlim_max > z2 : 
				zlim_max = z2
				dlim_max = d2
			if zlim_min < z1:
				zlim_min = z1
				dlim_min = d1

			print "Case 3"
		#	print "max,min", zlim_max, zlim_min, dlim_max,dlim_min
			for r in np.linspace(0,r_max,nz): # r in arcsec
				cr = sensitivity_function(r, mu=0, sigma=sigma)
				func_z = L_CO_prime_inv(l,co,cr,Sline,dv=200)
		
				#roots
				if  l < LCO_lim: 
					p = p1
					roots  = np.real((p-func_z).roots) #find solutions  and choose right root
					positive_roots =[ n for n in roots if n > 0]
				else  : 
					p = p2
					roots  = np.real((p-func_z).roots) #find solutions  and choose right root
					positive_roots =[ n for n in roots if n > 1.5]
				#print np.log10(l*cr),min(positive_roots)

				zlim=min(positive_roots)
				
				#print zlim
				if zlim > z2:
					dlim = d2
					zlim = z2
				elif zlim < z2 and zlim > z1:
					dlim = cosmo.luminosity_distance(zlim) #lower limit, defined by sensitivity
				elif zlim < z1: dlim = 0*u.Mpc

				#d1_prime =  (dlim_max-dlim_min)/nz * i #lower limit for a cilinder
				#print dlim
				h  = dlim - dlim_min # height of the element cilinder
				R = (np.radians(r/3600.) * dlim/ (1.+zlim)**2 )
				dR = (np.radians(dr/3600.) * dlim/ (1.+zlim)**2 )
				Volume += 2*np.pi*dR*R*h

				#print R,zlim,h
				#print "R,h,dR", R, h, dR,2*np.pi*dR*R*h
				#print dlim,h
				#else: 
				#	Volume += 0.


				i+=1
#choosing the regime based on the LCO' value. np.log10(LCO') < lim p1 > lim p2

print "zlim",z1,z2
print "LCOlim",  CO_Lum_prime2(z1,co, Sline, dv=200)/1.e8,CO_Lum_prime2(z2,co, Sline, dv=200)/1.e8
print "LCO",LCO/1e8
Z=np.linspace(0,10,1000)
co=CO[3]
Llim = np.log10(CO_Lum_prime2(1.5,co, Sline, dv=200)) # Luminocity limit for using different aproximaitons
print "Llim", Llim

plt.figure(1)
plt.axhline(z1)
plt.axhline(z2)



#looks like Llum > 1e10 has for J> 3 is always full volume size
for l in [1e6,1e7,1e8,1e9,1e10,1e11,1e12]:
	if np.log10(l) < Llim : 
		f = L_CO_prime_inv(l,co,1,Sline,dv=200) #get the value of dL^2/(1+z)
		roots  = np.real((p1-f).roots) #find solutions 
		positive_roots =[ n for n in roots if n > 0]
		print np.log10(l),positive_roots
		print min(positive_roots)
		zlim=min(positive_roots)
	
	else:
		f = L_CO_prime_inv(l,co,1,Sline,dv=200) #get the value of dL^2/(1+z)
		roots  = np.real((p2-f).roots) 
		positive_roots =[ n for n in roots if n > 1.5]
		print np.log10(l),min(positive_roots)
		zlim=min(positive_roots)
	plt.plot(np.log10(l),zlim, 'bo')
plt.show()

## caliculate the full cube volume (limits z1,z2 from frequency coverage)


for r in np.linspace(0,r_max,nz): # r in arcsec
	cr = sensitivity_function(r, mu=0, sigma=sigma)
	if  z2 < 1.5: 
		p = p1
	else:
		p = p2
	func_z = L_CO_prime_inv(LCO,co,cr,Sline,dv=200)
	
	root =  np.absolute((p-func_z).roots)
	#print root
	for s in root:
		
		if s  > z1 and s < z2 :
			zlim = s


	#print "zlim", zlim
	if zlim != 0:
		dlim = cosmo.luminosity_distance(zlim).value
		d1_prime = ddd * i #lower limit

		#plt.plot(r,d1_prime,'ro')
		#plt.plot(r,dlim,'bo')
		if d1_prime < d1:
			d1_prime = d1
		
		h  = dlim - d1_prime # height of the element cilinder

		zlim = z_at_value(cosmo.luminosity_distance, dlim*u.Mpc)
		if h > 0:
			R = (np.radians(r/3600.) * dlim/ (1.+zlim)**2 )
			dR = (np.radians(dr/3600.) * dlim/ (1.+zlim)**2 )
			Vol_max += 2*np.pi*dR*R*h

	i+=1
print Vol_max
#plt.axhline(d1)
#plt.show()

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

