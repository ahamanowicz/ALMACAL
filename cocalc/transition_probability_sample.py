import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
import astropy.units as u
from astropy import constants as const

mpl.rc ('xtick',labelsize=18)
mpl.rc ('ytick',labelsize=18)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#################################################

def CO_Lum_prime(z,fobs, Int_flux):
	
	dL = cosmo.luminosity_distance(z).value
	Lline_prime = 3.257e7 * Int_flux *  dL**2/(1+z)**3 / fobs**2 #K km/s pc^2

	return Lline_prime


def Mmol(Lco10,alpha=4.):
	#SLED choice 0- MW (normal SF, 1-M82 (ULIRG), 2 - Lagos+12 z=2 (SIMS) 3 - Lagos+12 z=2
	#returns LCO(1-0)' (K km/s pc^2) and  the H2 mas in Msol

	Mmol = alpha * Lco10
	return Mmol

def CO10_Lum_prime(z,Int_flux):
	co = 115.27

	dL = cosmo.luminosity_distance(z).value
	Lline_prime = 3.257e7 * Int_flux *  dL**2/(1+z)/co**2 #K km/s pc^2

	return Lline_prime

#################################################

######
# Take input flux + cenrtral frequency of the b=observed line and calculate the probability for the certain CO transition from Shark SAMS
######
print "################################################################################ \n Returns probability of the observed line to be a certain CO transition. \n Probabilities are based on the population of galaxies generated from Shark SAMS (Lagos et al 2018)\n for a galaxies with log(Molecular Mass) > 9.5.\n################################################################################"


## CO transitions
CO=[115.27,230.538,345.796,461.041,576.268,691.473,806.652, 921.7997,1036.912393,1151.985452] #GHz
CO_label = ["CO(1-0)","CO(2-1)","CO(3-2)","CO(4-3)","CO(5-4)","CO(6-5)","CO(7-6)", "CO(8-7)", "CO(9-8)", "CO(10-9)"]

Volume = 200000. #Mpc3 LAMACAL deep volume


#J "CO(1-0) - 0 (10) ","CO(2-1) - 1 (11) ","CO(3-2) - 2 (12)","CO(4-3) - 3 (13) ","CO(5-4) -4","CO(6-5) -5 ","CO(7-6) -6", "CO(8-7) -7", "CO(9-8) -8", "CO(10-9) -9"
# in Shark file 10, 11, 12 
cubes= np.loadtxt('candidates_list.txt', usecols=(1,),dtype='str',comments='#', unpack=True)
X,Y,Z,Width, SN, w50, f0, int_flux= np.loadtxt('candidates_list.txt', usecols=(3,4,5,8,6,10,9,11),dtype='float', unpack=True)
####
#print "central frequency [GHz]", freq0, "Integrated flux [mJy km/s]", flux
#### Create the table of possible transitions 
transition_prob= np.zeros((np.size(CO_label), 4))
#col 0 - J, col1 - z col2 - flux col3 - probability for the transition

##load transition bins
xedges= np.load("xbins.npy")
yedges=np.load("ybins.npy")
nbins=np.size(xedges)

## load completness

compeltnesH =np.load("completnes.npy")
completnes_x = np.load("compeltnes_xedges.npy") #SN
completnes_y = np.load("compeltnes_yedges.npy") #width km/s

Nprob,COprob,Z=[],[],[]
M=0
S = np.size(cubes)
transitions_array = np.zeros((S,10,5)) #[n_det][Jtrans][prob, L',z,L'(1-0),completnes]

#relaibility coeficient
rel_coeff = 0.9 #comes from sofia | to be improved

for s in range(S):

	flux = int_flux[s] #mJy
	freq0 = f0[s]
	
	Nprob,COprob,Z,LCO10=[],[],[],[]
	M=0

	#check completnes for detection parametes \
	sn, w  = SN[s], Width[s]
	#print sn,w
	if w < 50. : w= min(completnes_y )
	for k in range(np.size(completnes_x)-1):
		x1,x2 = completnes_x[k],completnes_x[k+1]

		if x1 <= sn <= x2: 
			x=k
			for l in range(np.size(completnes_y)-1):

				y1,y2 = completnes_y[l],completnes_y[l+1]
				#print y1,y2
				if y1 <= w <= y2: y=l
	
	comp_coef = compeltnesH[y][x] #completnes soefficient for a given detections 1/coeffishent is the actual nr of objects like that
	#print "compl",comp_coef

	for j in range(np.size(CO_label)):
		jup=j+1
		z = CO[j]/freq0 -1.
		flu = np.log10(flux)
		# read in the histogram 
		H = np.load(CO_label[j]+"_hist.npy")
		M += np.sum(H)
		if jup > 1 :
			LCO10_H = np.load("I_CO10_"+str(jup)+".npy")
		#print z, flu
		if z > 0 and z < 10:
			#print z
			#check which bin you are in for k in range(nbins-1):
			for k in range(nbins-1):
				#print xedges[k],xedges[k+1], yedges[l],yedges[l+1],H.T[k][l]
				x1,x2 = xedges[k],xedges[k+1]
				if x1<=z <=x2: 
					x = k
				
					for l in range(nbins-1):
						y1,y2 = yedges[l],yedges[l+1]
						if y1<=flu <=y2: 
							
							y = l
							
							Nprob.append(H[y][x])
							COprob.append(j)
							Z.append(z)
							if jup > 1:
								Ico10 = LCO10_H[y][x] #it is a log
								lco10_prime = CO10_Lum_prime(z,(10**Ico10)/1.e3) #from mJy to Jy
								LCO10.append(lco10_prime)
							else:
								
								lco10_prime = CO10_Lum_prime(z,flux/1.e3)
								LCO10.append(lco10_prime)
							#print H.T[y][x], j 
		elif z < 0:
			COprob.append(j+1)
			Z.append(str(0))
			Nprob.append(str(0))
			LCO10.append(str(0))
	Nprob = np.array(Nprob, dtype='float')
	N = np.sum(Nprob)	
	Mnew=np.sum(Nprob/M)
	Nprob_all = Nprob/M
	LCO10 = np.array(LCO10,dtype='float')
	#print Nprob, LCO10
	#f.write(cubes[s]+ " "+str(freq0)+" "+str(int_flux[s]))
	for i in range(np.size(COprob)):
		j = COprob[i]
		z = CO[j]/freq0 -1.
		prob = Nprob[i]/M/Mnew
		jup = int(j) +1
		if prob > 0: 
			Lprime = CO_Lum_prime(z,freq0, flux)
			
			L10prime = LCO10[i]
			print "transition", CO_label[j], "z =", np.round(z,3), "probability", np.round(prob*100,3), "log(Lprime)", np.round(np.log10(Lprime),3),"log(L10prime)", np.round(np.log10(L10prime),3) 
			#f.write("transition", CO_label[j], "z =", np.round(z,3), "probability", np.round(prob*100,3), "log(Lprime)", np.round(np.log10(Lprime),3), "log(Mmol)",np.round(np.log10(MH2),3),"log(L10prime)", np.round(np.log10(L10prime),3) )
		else: 
			prob , Lprime,L10prime= 0,0,0
		array=[prob,Lprime,z,L10prime,comp_coef]
		#results_array[l][n][:]
		#print "array", array
		for k in range(5):
			array[k]
			transitions_array[s][i][k] = array[k]


### use transition array to draw the probable samples
ntry= 1000 # muber of tries for probability drawing
results_array=np.zeros((ntry,S,3)) #[ntries],[nprob][z,L'1-0,MH2, comp]
nr_detections= np.shape(transitions_array)[0]
nr_trans = np.shape(transitions_array)[1]

####create redhsift bins for roH2 - 

## !!!! TEST !!!
# ASPECS - like bins egdes
z_bins = np.array([0,0.5,1.0,2.0,3.,4.])
# corresponding volume: in these redhsift bins sum the volumes coming from each transition to get Volume [Mpc^3] per z bin
Volumes  = np.array([20, 20.,20.,20.,20]) #bvolume per redshfi bin

# luminocity bins 0.5 dex
luminosity_bins = np.arange(8,12.5,0.5)

LFa1, LFa2, LFa3, LFa4, LFa5=np.zeros((1,np.size(luminosity_bins)-1)), np.zeros((1,np.size(luminosity_bins)-1)),np.zeros((1,np.size(luminosity_bins)-1)),np.zeros((1,np.size(luminosity_bins)-1)),np.zeros((1,np.size(luminosity_bins)-1))
LFall = np.zeros((1,np.size(luminosity_bins)-1))
for l in range(ntry): #realsations of the sample
	lco10_sample, zsmample=[],[]
	for n in range(nr_detections):
		P,prob = [],[]
		#for each detection
		for m in range(nr_trans):
			# 1) create a table of transitions probabilities and table of indexes
			
			P = np.append(P,transitions_array[n][m][0]) 
		Indexes = np.array(range(10))
		#print P
		# 2) draw a transition to ba assigned to the detection, according to the probabilities
		#draw a transition index I
		I = np.random.choice(Indexes,size=1,p=P)[0]

		#print I, transitions_array[n][I][:]
		for k in range(3):
			results_array[l][n][k] = transitions_array[n][I][k+2] #[z,L'(1-0), completnes]

	#general LF
	ncounts_all= np.zeros((1,np.size(luminosity_bins)-1))
	for j in range(np.size(luminosity_bins)-1):

		for i in range(nr_detections):
			lum = results_array[l][i][1]
			comp = results_array[l][i][2]
			#print lum,comp
			if lum > 0:
				#print np.log10(lum)
				if luminosity_bins[j] <= np.log10(lum) <  luminosity_bins[j+1]:
					
					ncounts_all.T[j] += 1./comp
				
	LFa = ncounts_all*0.9/Volume		
	LFall = np.vstack((LFall,LFa))	

	#LCO'(1-0) sample, divide into z bins
	# lco bins, get indexes of the detections from the bigger array
	l1,l2,l3,l4,l5=np.array([[0,0]]),np.array([[0,0]]),np.array([[0,0]]),np.array([[0,0]]),np.array([[0,0]])
	#divide ito redshift bins
	z1,z2  = z_bins[0],z_bins[1]
	for n in range(nr_detections):
			if  z_bins[0] <= results_array[l][n][0] < z_bins[1]:
				l1 = np.vstack((l1, [np.log10(results_array[l][n][1]),results_array[l][n][2]]))

			elif z_bins[1]  <= results_array[l][n][0] < z_bins[2]:
				l2 = np.vstack((l2, [np.log10(results_array[l][n][1]),results_array[l][n][2]]))
			elif z_bins[2]  <= results_array[l][n][0] < z_bins[3]:
				l3 = np.vstack((l3, [np.log10(results_array[l][n][1]),results_array[l][n][2]]))			
			elif z_bins[3]  <= results_array[l][n][0] < z_bins[4]:
				l4 = np.vstack((l4, [np.log10(results_array[l][n][1]),results_array[l][n][2]]))			
			elif z_bins[4]  <= results_array[l][n][0] < z_bins[5]:
				l5 = np.vstack((l5, [np.log10(results_array[l][n][1]),results_array[l][n][2]]))			
#	print np.shape(l1)[0]

	# CO(1-0) LF for each z bin 
	# luminoity finction for each redshfit bin, luminocity bins every 0.5 dex

	# bin 1 0- 0.5
	#if np.shape(l1)[0] > 1:
	#	for j in np.size(luminosity_bins-1):
	#		print 0
	#bin 0.5 - 1

	vol = 20. ## !!!!!! repair	
	if np.shape(l1)[0] > 1:
		ncounts= np.zeros((1,np.size(luminosity_bins)-1))
		for j in range(np.size(luminosity_bins)-1):

			for i in range(np.size(l2.T[0])):
				if l1.T[0][i] > 0:
					if luminosity_bins[j] <= l1.T[0][i] <  luminosity_bins[j+1]:
						
						ncounts.T[j] += 1./l1.T[1][i]
		LF1 = ncounts*0.9/vol		
		LFa1 = np.vstack((LFa1,LF1))	
	if np.shape(l2)[0] > 1:
		ncounts = np.zeros((1,np.size(luminosity_bins)-1))
		for j in range(np.size(luminosity_bins)-1):

			for i in range(np.size(l2.T[0])):
				if l2.T[0][i] > 0:
					if luminosity_bins[j] <= l2.T[0][i] <  luminosity_bins[j+1]:
						
						ncounts.T[j] += 1./l2.T[1][i]
	
		LF2 = ncounts*0.9/vol	
		LFa2 = np.vstack((LFa2,LF2))	
	if np.shape(l3)[0] > 1:
		ncounts = np.zeros((1,np.size(luminosity_bins)-1))
		for j in range(np.size(luminosity_bins)-1):

			for i in range(np.size(l3.T[0])):
				if l3.T[0][i] > 0:
					if luminosity_bins[j] <= l3.T[0][i] <  luminosity_bins[j+1]:
						
						ncounts.T[j] += 1./l3.T[1][i]
	
		LF3 = ncounts*0.9/vol
		LFa3 = np.vstack((LFa3,LF3))	
	if np.shape(l4)[0] > 1:
		ncounts = np.zeros((1,np.size(luminosity_bins)-1))
		for j in range(np.size(luminosity_bins)-1):

			for i in range(np.size(l4.T[0])):
				if l4.T[0][i] > 0:
					if luminosity_bins[j] <= l4.T[0][i] <  luminosity_bins[j+1]:
						ncounts.T[j] += 1./l4.T[1][i]
	
		LF4 = ncounts*0.9/vol
		LFa4 = np.vstack((LFa4,LF4))	
	if np.shape(l5)[0] > 1:
		ncounts = np.zeros((1,np.size(luminosity_bins)-1))
		for j in range(np.size(luminosity_bins)-1):

			for i in range(np.size(l5.T[0])):
				if l5.T[0][i] > 0:
					if luminosity_bins[j] <= l5.T[0][i] <  luminosity_bins[j+1]:
						
						ncounts.T[j] += 1./l5.T[1][i]
	
		LF5 = ncounts*0.9/vol
		LFa5 = np.vstack((LFa5,LF5))	

### mean and sigmas of LF 
plt.figure(1)

for i in range(np.size(LFall[0])):
	#print LFa2.T
	#print np.mean(LFall.T[i]), np.std(LFall.T[i])

	if np.mean(LFall.T[i]) > 0 :
		plt.errorbar(luminosity_bins[i], np.log10(np.mean(LFall.T[i])), fmt='o',c='navy')

	
	plt.plot(luminosity_bins[i], np.log10(np.max(LFall.T[i])), 'ro')
	plt.plot(luminosity_bins[i], np.log10(np.min(LFall.T[i])), 'ro')
	plt.ylim([-6,1])
	plt.xlim([7.5,12.5])
plt.show()
"""
plt.figure(1, figsize=(8,8))
Mh2_array=[results_array.T[:][0][l],results_array.T[:][2][l]]
for l in range(ntry):
	array = np.array([results_array.T[:][0][l],results_array.T[:][2][l]])
	array_sort = np.sort(array,axis=1)
	print array_sort[0]
	plt.plot(array_sort[0],np.log10(array_sort[1]/Volume))
	Mh2_array=np.append(Mh2_array,array,axis=1)

Mh2_array_sort = np.sort(Mh2_array)
print Mh2_array_sort[0]

plt.show()
"""
	