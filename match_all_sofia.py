import matplotlib.pyplot as plt
import numpy as np
from scipy import constants as const
from scipy.stats import skewnorm
from scipy.integrate import quad
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle, Ellipse
from scipy.signal import gaussian, fftconvolve
from astropy.io import fits
from astropy.wcs import WCS
import astropy
from scipy.ndimage.interpolation import rotate
import scipy.integrate as integrate
import sys
import os


import matplotlib as mpl

mpl.rc ('xtick',labelsize=15)
mpl.rc ('ytick',labelsize=15)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

red='#d11141'
yellow='#ffc425'
blue='#00aedb'
green='#00b159'
orange='#f37735'
#fname = sys.argv[1]
#cubes= np.loadtxt(fname, usecols=(0,),dtype='str')
matches_array=[]
Width, SN = np.loadtxt('candidates_list.txt', usecols=(8,6),dtype='float', unpack=True)
clas= np.loadtxt('candidates_list.txt', usecols=(12,),dtype='str',comments='#', unpack=True)

fig = plt.figure(1,figsize=(8,8))

sn_a,w_a = [],[]
sn_m, w_m=[],[]

cubes,flags = np.loadtxt("fits_list_all.in", dtype='str', usecols=(0,1), delimiter=' ', unpack=True)
N_cubes = np.size(cubes)

for i in range(N_cubes):
	pref = cubes[i]
	flag = flags[i]
	#print flag, map(int,flag.split(','))
	if flag != '':
		flag = map(int,flag.split(','))
	else: flag =[]
	#print pref
	data=np.loadtxt(pref.replace(".fits","_mock_nozero_cat.ascii"),comments='#',dtype='str')
	mock = np.loadtxt(pref.replace(".fits","_mock.txt"),skiprows=2, dtype='str')
	mock_s = np.loadtxt(pref.replace(".fits","_mock.txt"), usecols=(0,6,7),skiprows=2)
	
	M =  np.shape(data)[0]
	N = np.shape(mock)[0]

	with open(pref.replace(".fits", ".stats")) as f:
		content = f.read().splitlines()
	a,b = content[3].split('\t')[1].split(' ')
	a,b =float(a),float(b)
	Mock,Match, Match_mock=[],[],[]
	f1,f2 = float(content[7].split()[2]), float(content[7].split()[3])
	Zmax = float(content[6].split()[1])
	f0 = (f1+f2)/2.
	dZmax = 700*1.e3/const.c *f0 #in GHz
	dZmax = dZmax * Zmax / abs(f2-f1)
	#print  dZmax, Zmax
	g = open(pref.replace("fits","matches"),'w')
	for n in range(N):
		Z_width, Zm = int(mock[n][8]),int(mock[n][3])
		#check if mock not flagged
		z1,z2= int(Zm - Z_width/2.), int(Zm + Z_width/2.)
		mock_z =  range(z1,z2,1)
		overlap =  [x for x in mock_z if x in flag]
		if 1==1: #not (len(mock_z) > 5. and len(overlap) > 2) or len(overlap) < 0.:
			sn_m.append(mock[n][6])
			w_m.append(mock[n][7])

	for m in range(M):
		line=''
		#if data[m][28] >= 3.:
		sn = round(float(data[m][26])/float(data[m][30]),3) # f_peak / rms THIS IS WORNG STILL
		f_peak = float(data[m][26]) *1.e3 #mJy/beam
		X,Y,Z = float(data[m][5]),float(data[m][6]),float(data[m][7])
		N_chan = int(data[m][15])
	#	line=str(m+1) +' X '+ data[m][5] + ' Y ' + data[m][6] + ' Z ' + data[m][7] +' S/N '+str(sn)+' width '+str(data[m][15])
		for n in range(N):
			Xm,Ym,Zm = int(mock[n][1]),int(mock[n][2]),int(mock[n][3])
			Z_width = int(mock[n][8])
			f_peak_mock, width_kms_mock = float(mock[n][5]),float(mock[n][7])
			#check if mock not flagged
			z1,z2= int(Zm - Z_width/2.), int(Zm + Z_width/2.)
			mock_z =  range(z1,z2,1)
			
			overlap =  [x for x in mock_z if x in flag]
			
			if 1==1: #not (len(mock_z) > 5. and len(overlap) > 2) or len(overlap) < 0.:
				
				# matching Z and sdistance
				dist = np.linalg.norm(np.array([X,Y]) - np.array([Xm,Ym]))
				if dist <= np.sqrt(a**2 + b**2): # np.sqrt(3^2 + 3^2)
					Z_dist = abs(Zm-Z)
					if Z_dist <= dZmax:
						
						if str(n) not in Mock:
							#print n
							Mock.append(str(n))

							#plt.plot(mock[n][6],mock[n][7],'ro',alpha=0.5)
						if str(n) not in Match_mock:
							Match_mock.append(str(n))
							sn_a.append(mock[n][6])
							w_a.append(mock[n][7])
						Match.append(str(m))
						line=str(m+1) +' X '+ data[m][5] + ' Y ' + data[m][6] + ' Z ' + data[m][7] +' S/N '+str(sn)+' width '+str(data[m][15])+' f_peak ' + str(f_peak)+ ' match '+str(n+1)+ ' S/N ' + str(mock[n][6]) + ' width ' +str(mock[n][8])+ ' f_peak ' + str(f_peak_mock)+"\n"
						plt.scatter(float(mock[n][6]),float(mock[n][7]),edgecolor=red,facecolor='none',linewidth=6)
					
			
			#line+='\n'	

		#if data[m][31] > 3.:
		g.write(line)	
	#plt.legend()

	g.close()

	#print candidates
	h = open(pref.replace('fits','cand'),'w')
	for m in range(M):
		sn = round(float(data[m][26])/float(data[m][30]),3)
		line=''
		if str(m) not in Match:

			Z = data[m][7]
			line=str(m+1) +' X '+ data[m][5] + ' Y ' + data[m][6] + ' Z ' + data[m][7] +' S/N '+str(sn)+' width '+str(data[m][15]) + '\n'
			
			if float(Z) > 3. and float(Z) < Zmax -3.:
				
				h.write(line)
	h.close()
	matches = np.size(Mock)
	print pref, "detected", M,"matched", matches, "multiple matches", np.size(Match)-matches

	#print matches_array
plt.plot(sn_m,	w_m,'o',color='navy',ms=5)
#add sources missed form matching (reasons, being tmaerged detections of too closely injected sources, source injected at Z=0 and being shifted etc)
sn_missed, w_missed = np.loadtxt('missing_mathced_sources.txt', skiprows=2, dtype='float', delimiter=',', usecols=(8,9), unpack=True)
plt.scatter(sn_missed,w_missed,edgecolor=red,facecolor='none',linewidth=6)
plt.plot(sn_missed,w_missed,'o',color='navy',ms=5)
plt.axvline(3,lw=3,c=yellow)
plt.xlim([0.8, 8.2])
plt.ylim([40,820])
plt.grid("on")
plt.xlabel("S/N", fontsize=20)
plt.ylabel("width [km/s]", fontsize=20)

fig2 = plt.figure(2,figsize=(10,10))
fig2.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.23)
sn_a = np.array(sn_a, dtype='float')
sn_a = np.append(sn_a, sn_missed,axis=0)
sn_m = np.array(sn_m, dtype='float')
sn_m = np.append(sn_m, sn_missed, axis=0)
w_a = np.array(w_a, dtype='float')
w_a = np.append(w_a, w_missed, axis=0)
w_m = np.array(w_m, dtype='float')
w_m = np.append(w_m, w_missed, axis=0)


ax1 = plt.subplot2grid((3,3), (0,1), rowspan=2, colspan=2)

H,xedges,yedges = np.histogram2d(sn_a,w_a,bins=[10,10])
H_all,xedges_all,yedges_all = np.histogram2d(sn_m,w_m,bins=[xedges,yedges])
H_frac = H/H_all
print H_frac
grid = H_frac.T

X,Y=np.meshgrid(xedges,yedges)
im=plt.pcolormesh(X,Y,grid)

for i in range(np.size(clas)):	
	if clas[i] == "YES": c = "red"
	else : c = 'yellow'
	plt.scatter(SN[i], Width[i], c=c, s=50)
N,M = range(np.size(w_a)), np.size(sn_a)
#plt.axis("off")
#plt.colorbar(orientation='horizontal')
print M
cb_ax = fig2.add_axes([0.83, 0.38, 0.02, 0.52])
cbar = fig2.colorbar(im, cax=cb_ax)
cbar.set_clim(0,1.0)
ax2 = plt.subplot2grid((3,3),(2,1), colspan=2)

hist_sn_m, bins_sn_m = np.histogram(sn_m, bins=10)
hist_sn_det, bins_sn_det = np.histogram(sn_a, bins=bins_sn_m)
bin_width = bins_sn_m[1] - bins_sn_m[0]

#all SN all widths colalpesed
fraction=np.array(hist_sn_det,dtype='float')/np.array(hist_sn_m, dtype='float')
plt.plot(bins_sn_m[1:]- bin_width/2.,fraction, 'o--', c='navy')
#fractions of the 2d histogram 
## width > 400

plt.ylim([0.0,1.0])
plt.xlim([xedges[0], xedges[-1]])
plt.xlabel("S/N", fontsize=20)
plt.ylabel("Detection fraction", fontsize=20)
ax3 = plt.subplot2grid((3,3),(0,0), rowspan=2)

hist_w_m, bins_w_m = np.histogram(w_m, bins=10)
hist_w_det, bins_w_det = np.histogram(w_a, bins=bins_w_m)
bin_width = bins_w_m[1] - bins_w_m[0]
fraction=np.array(hist_w_det,dtype='float')/np.array(hist_w_m, dtype='float')
plt.plot(fraction,bins_w_m[1:]-bin_width/2., 'o--', c='navy')
plt.xlim([1.0,0.0])
plt.ylabel("width [km/s]", fontsize=20)
plt.xlabel("Detection fraction", fontsize=20)
plt.ylim([yedges[0], yedges[-1]])

#plt.tight_layout()
fig.savefig("Mock_detected_all_candidates.pdf")
fig2.savefig("Detection_fraction.pdf")
plt.show()
