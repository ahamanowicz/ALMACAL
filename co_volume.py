import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
import astropy.units as u
from astropy import constants as const
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
	Lline_prime = 3.47e7 * Sline * dv *  dL**2/(1+z)**3 / f_obs**2
	return Lline_prime

def L_CO_prime_inv(LCO,fCO,cr,Sline,dv=200):
	func = LCO*cr*fCO**2/(3.47e7 * Sline * dv )
 	return func

def sensitivity_function(r, mu, sigma):

	I = np.exp(-np.power(r - mu, 2.) / (2 * np.power(sigma, 2.)))

	return I

def E(z):
	#Planck15 cosmology
	Om_M = cosmo.Om0 #Omega matter
	Om_L = 1-Om_M # omega Lambda , since Omega_k = 0
	#Om_k = 1- Om_M - Om_L # = 0 flat Universe

	E = np.sqrt(Om_M*(1+z)**3 + Om_L)
	return E

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
# 0-CO(1-0) 1-CO (2-1) 2-CO(3-2) 3-CO(4-3) 4-CO(5-4) 5-CO(6-5) 6-CO(7-6) 7-CO(8-7) 8-CO(9-8)

CO = np.array([115.27,230.538, 354.796, 461.041, 576.268, 691.473, 806.652,921.7997,1036.9124])

# transition index 
k = 1 # CO(3-2)

####### read data from the cube *stats files
filename='good_cubes.txt'#'volumes.in'
files = np.loadtxt(filename, usecols=(0,),dtype='str')
#cosmology
c =  const.c.to('km/s')
H0 =  cosmo.H(0)
DH = c/H0
#luminocity range
lum_n = 50
N = np.size(files)
logLUM = np.linspace(5,13,lum_n)#luminocity sampling
LUM = 10**logLUM
M = np.size(logLUM)

plt.figure(1,figsize=(13,10))
g= open('Volumes_CO.txt', 'w')
g.write("# CO_trans number, zmean, Vol_total\n")
g.write("# 0-CO(1-0) 1-CO (2-1) 2-CO(3-2) 3-CO(4-3) 4-CO(5-4) 5-CO(6-5) 6-CO(7-6) 7-CO(8-7) 8-CO(9-8)\n")
g.write("# Vol_total for logLUM np.linespace(5,13,50)\n")

for k in range(np.size(CO)):

	L_VOL = np.zeros((lum_n,N+2))
	L_VOL.T[0] = logLUM

	Vol_all_cubes=np.zeros((lum_n))
	Vol_all_cubes_z005=np.zeros((lum_n))
	Vol_all_cubes_z051=np.zeros((lum_n))
	Vol_all_cubes_z115=np.zeros((lum_n))
	Vol_all_cubes_z152=np.zeros((lum_n))
	Vol_all_cubes_z2 = np.zeros((lum_n))
	for i in range(0,N):
		L_VOL = np.zeros((lum_n,N+2))
		L_VOL.T[0] = logLUM
		fname = files[i]
		with open("cubes_stats/"+fname) as f:
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
		r_max = np.sqrt(2)*sigma #arcsec | sqrt(2) ~ 1.7
		#print imsize_radians,sigma
		Sline = 5 *rms # line strength
		dv = 200 # [km/s] width of the line 

		##### redshift coverage of the transition  [cube limits] #######

		co = CO[k] # CO transition central frequency [GHz]
		z1,z2 = co/f2-1, co/f1 -1 #z1 lower, z2 higher redshift | limiting redshft from the band coverage

		# security 
		if z2 > 0:
			print fname
			print z1,z2
			if z1 < 0 : z1=0
			d1 = cosmo.luminosity_distance(z1) # lower distance limit [Mpc]
			d2 = cosmo.luminosity_distance(z2) #higher distance limit [Mpc]
			#### INTEGRATION ####
			#all cosmological parameters in LUMINOCITY 
			nz = 100.
			dd = (z2-z1)/nz
			#print "r_max",r_max, "z1,z2",z1,z2
			dr = r_max/nz

			## maximum value of the volume - turnicated cone -> in the end I am using the cyllinder
			#for the R estimation using mean redshift
			zmean = (z2+z1)/2.
			p = cosmo.kpc_proper_per_arcmin(zmean) #conversion factor kpc/arcmin
			#angular size to physical
			R = (p*(r_max/60.)*u.arcmin).to(u.Mpc) #Mpc 
			Rc = R*(1+zmean) #to comoving 
			print R,Rc,d1,d2
			#Vol_max = 1./3. *np.pi * abs(d2-d1) *(R2**2 + R2*R1 + R1**2) #cone
			Vol_max_2 = np.pi * abs(d2-d1) * R**2 #cylinder, based on luminocity, not comoving
			#print "R1,R2,d1,d2,z1,z2",R1,R2,d1,d2, z1,z2
			print "Volume max", Vol_max_2
			
			#comoving volume correction
			Om_beam = np.pi * np.radians(r_max/3600.)**2
			
			element = 0.
			nz = 100.
			dz = (z2-z1)/nz

			for z in np.linspace(z1,z2,nz):
				dL = cosmo.luminosity_distance(z)

				element+= dL**2 / (1+z)**2 / E(z) * dz

			Vc_int = element*DH * Om_beam

			alpha = Vc_int/Vol_max_2
			
			#chosing the regime

			LCO_lim = CO_Lum_prime2(1.5,co, Sline, dv=200)
			#print "LCOlim", np.log10(LCO_lim), z1,z2
			#print d1,d2

			for j in range(0,M):

				l = LUM[j]
				Volume = 0.
				
				if  l < LCO_lim: 
					p = p1
					#print "p1"
					###limiting redshifts from sensitivity
					#max
					cr=1
					func_z = L_CO_prime_inv(l,co,cr,Sline,dv=200)
					roots  = np.real((p-func_z).roots) #find solutions  and choose right root
					positive_roots =[ n for n in roots if n > 0]
					zlim_max = min(positive_roots) #in the center of an image
					#print "max",positive_roots
					#min
					cr = sensitivity_function(r_max, mu=0, sigma=sigma)
					func_z = L_CO_prime_inv(l,co,cr,Sline,dv=200)
					roots  = np.real((p-func_z).roots) #find solutions  and choose right root
					positive_roots =[ n for n in roots if n > 0]
					zlim_min = min(positive_roots) #at the edge of he image
					#print "min",positive_roots
					if (zlim_min > zlim_max) : zlim_min,zlim_max = zlim_max,zlim_min
					dlim_max = cosmo.luminosity_distance(zlim_max)
					dlim_min = cosmo.luminosity_distance(zlim_min)
					#print "zmax, zmin, LCO", zlim_max,zlim_min, np.log10(l)

				else: 
					p = p2
					#print "p2"
					cr=1
					func_z = L_CO_prime_inv(l,co,cr,Sline,dv=200)
					roots  = np.real((p-func_z).roots) #find solutions  and choose right root
					positive_roots =[ n for n in roots if n > 1.5]
					zlim_max = min(positive_roots)
					#print "max", positive_roots
					#min
					cr = sensitivity_function(r_max, mu=0, sigma=sigma)
					func_z = L_CO_prime_inv(l,co,cr,Sline,dv=200)
					roots  = np.real((p-func_z).roots) #find solutions  and choose right root
					positive_roots =[ n for n in roots if n > 1.5]
					zlim_min = min(positive_roots)
					#print "min", positive_roots
					if (zlim_min > zlim_max) : zlim_min,zlim_max = zlim_max,zlim_min
					dlim_max = cosmo.luminosity_distance(zlim_max)
					dlim_min = cosmo.luminosity_distance(zlim_min)
				
				#print "zmax, zmin, LCO", zlim_max,zlim_min, np.log10(l)
				#print "dmax, dmin", dlim_max,dlim_min
				#print "z1,z2", z1,z2
				
				print zlim_max,zlim_min, z1,z2
				### CASE 1 - zlim_max > z2 and zlim_min  > z2; volume -> volume of the truncated cone
				if zlim_max > z2 and zlim_min  > z2:
					print "Case 1",np.log10(l)
					Volume = Vol_max_2

				### CASE 2 - zlim_max < z1
				elif zlim_max < z1 :
					print "Case 2",np.log10(l)
					Volume = 0*(u.Mpc)**3
				### CASE 3 - zlim_max > z2 and zlim_min < z2  (always check if zlim_min > z1)

				elif zlim_max > z2 and zlim_min < z2:	

					print "Case 3",np.log10(l)
					#print "z1,z2,zlimax,zlimmin",z1,z2,zlim_max,zlim_min
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
							p = cosmo.kpc_proper_per_arcmin(z2)
							R = (p*(r/60.)*u.arcmin).to(u.Mpc)
							dR = (p*(dr/60.)*u.arcmin).to(u.Mpc)
							#R = (np.radians(r/3600.) * d2/ (1.+z2)**2 )
							#dR = (np.radians(dr/3600.) * d2/ (1.+z2)**2 )
							Volume += 2*np.pi*dR*R*h	
							#print d2,d1,dlim,Volume
						elif zlim > z1:
						
							h = dlim - d1
							#print d2,d1,dlim
							p = cosmo.kpc_proper_per_arcmin(zlim)
							R = (p*(r/60.)*u.arcmin).to(u.Mpc)
							dR = (p*(dr/60.)*u.arcmin).to(u.Mpc)
							Volume += 2*np.pi*dR*R*h
						else:
							Volume +=0*(u.Mpc)**3
						#print zlim,z2,z1, Volume
					#print Volume*alpha		

				### CASE 4 - zlim_max < z2  and zlim_min > z1  (always check if zlim_min > z1)
				elif zlim_max < z2 and zlim_min  > z1:
					print "Case 4",np.log10(l)
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
						p = cosmo.kpc_proper_per_arcmin(zlim)
						R = (p*(r/60.)*u.arcmin).to(u.Mpc)
						dR = (p*(dr/60.)*u.arcmin).to(u.Mpc)
						#R = (np.radians(r/3600.) * dlim/ (1.+zlim)**2 )
						#dR = (np.radians(dr/3600.) * dlim/ (1.+zlim)**2 )
						Volume += 2*np.pi*dR*R*h
						
						#Volume -=1./3.*np.pi*(np.radians(r_max/3600.) * dlim_min/ (1.+zlim_min))**2*dlim_min 
				### CASE 5 - zlim_max < z2  and zlim_min < z1  (always check if zlim_min > z1)
				elif  zlim_max < z2 and zlim_max > z1 and zlim_min  < z1:
					print "Case 5",np.log10(l)
					
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
						#print positive_roots
						#print np.log10(l*cr),min(positive_roots)
						zlim=min(positive_roots)
						dlim = cosmo.luminosity_distance(zlim)
						
						if zlim > z1:
							
							#d1_prime = dlim_max/nz * i #lower limit
							h  = dlim - d1 # height of the element cilinder
							p = cosmo.kpc_proper_per_arcmin(zlim)
							R = (p*(r/60.)*u.arcmin).to(u.Mpc)
							dR = (p*(dr/60.)*u.arcmin).to(u.Mpc)
							#R = (np.radians(r/3600.) * dlim/ (1.+zlim)**2 )
							#dR = (np.radians(dr/3600.) * dlim/ (1.+zlim)**2 )
							Volume += 2*np.pi*dR*R*h
							#print dlim,d1,h
							#print zlim,z1
						else: 
							Volume += 0.*(u.Mpc)**3
					#	print zlim,Volume
				
				#print Volume
				Volume_c = Volume * alpha #convertion to comoving 
				Volume_c = Volume_c.value
				#newrow=[np.log10(l),Volume_c] #correction to comoving
				#L_VOL.append(newrow)
				#print Volume_c
				L_VOL.T[i+1][j] = Volume_c

			
			Vol_total = np.sum(L_VOL.T[1:N+1], axis=0)
			Vol = np.array(Vol_total,dtype='str')
			Vol = " ".join(Vol)
			g.write(str(k)+' '+str(zmean)+ ' ')
			g.write(Vol)
			g.write('\n')
		 	#print zmean
			print Vol_total,zmean
			if zmean <= 0.5:
				Vol_all_cubes_z005+=Vol_total
			elif zmean > 0.5 and zmean <=1.0:
				Vol_all_cubes_z051+=Vol_total
			elif zmean > 1.0 and zmean <= 1.5:
				Vol_all_cubes_z115+=Vol_total
			elif zmean >1.5 and zmean <=2.0:
				Vol_all_cubes_z152+=Vol_total
			else:
				Vol_all_cubes_z2 +=Vol_total
			Vol_all_cubes+=Vol_total

	
		
	plt.subplot(3,3,k+1)
	plt.plot(logLUM, Vol_all_cubes_z005, '-ro', label='z=[0-0.5]')#,c='b')
	plt.plot(logLUM, Vol_all_cubes_z051, '-go', label='z=[0.5-1.0]')#,c='b')
	plt.plot(logLUM, Vol_all_cubes_z115, '-bo', label='z=[1.0-1.5]')#,c='b
	plt.plot(logLUM, Vol_all_cubes_z152, '-mo', label='z=[1.5-2.0]')#,c='b')
	plt.plot(logLUM, Vol_all_cubes_z2, '-ko', label='z > 2 ')#,c='b')
	plt.legend()
	co_trans = 'CO('+str(k+1)+'-'+str(k)+')'

	plt.annotate(co_trans, xy=(0.05, 0.25), xycoords='axes fraction',fontsize=15)
	plt.xlabel("log10(LCO_prime)")
	plt.ylabel("Volume comoving [Mpc^3]")
g.close()
plt.savefig("Volumes.pdf")
plt.show()
"""	
#print Vol_all_cubes
plt.figure(2,figsize=(10,8))
plt.plot(logLUM, Vol_all_cubes, 'go')
#plt.plot(logLUM, Vol_all_cubes, '-yo')

#g= open(filename.replace('.txt','_CO'+str(k+1)+str(k)+'.vol'), 'w') 
#g= open('_CO'+str(k+1)+str(k)+'.vol', 'w') 
#g.write("logCO_LUM, Volume [Mpc^3]\n")
#for h in range(np.size(logLUM)):
#	g.write(str(float(round(L_VOL.T[0][h],6)))+' '+str(float(round(Vol_total[h],6)))+'\n')
#
"""
