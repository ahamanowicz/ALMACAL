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


cubes= np.loadtxt('fits_list.in', usecols=(0,),dtype='str')
matches_array=[]
fig = plt.figure(1,figsize=(8,8))
sn_a,w_a = [],[]

Zmax = 700 #km/s


for pref in cubes:
	#print pref
	data=np.loadtxt(pref.replace(".fits","_results.txt"),comments='#',dtype='str')

	mock = np.loadtxt(pref.replace(".fits","_mock.txt"),skiprows=2, dtype='str')
	mock_s = np.loadtxt(pref.replace(".fits","_mock.txt"), usecols=(0,6,7),skiprows=2)
	#print mock_s.T[1]
	plt.plot(mock_s.T[1],mock_s.T[2],'bo')
	M =  np.shape(data)[0]
	N = np.shape(mock)[0]
	for i in np.arange(0,9,1):
		plt.axvline(i,c='grey',linewidth=1)
	for i in range(0,900,100):
		plt.axhline(i,c='grey',linewidth=1)

	stats = pref.replace("fits","stats")
	with open(stats) as f:
		content = f.read().splitlines()
	a,b = content[3].split('\t')[1].split(' ')
	a,b =float(a),float(b)
	Mock,Match=[],[]

	g = open(pref+".out",'w')
	h = open(pref+'.cand','w')
	for m in range(M):
		line=''
		#if data[m][28] >= 3.:
		sn = data[m][20]
		
		X,Y,Z = float(data[m][31]),float(data[m][32]),float(data[m][33])
		N_chan = int(data[m][28])
		line=str(m+1) +' X '+ data[m][31] + ' Y ' + data[m][32] + ' Z ' + data[m][33] +' S/N '+str(sn)+' width '+str(data[m][28])
		for n in range(N):
			Xm,Ym,Zm = int(mock[n][1]),int(mock[n][2]),int(mock[n][3])
			Z_width = int(mock[n][8])
			dist = np.linalg.norm(np.array([X,Y]) - np.array([Xm,Ym]))
			if dist <= np.sqrt(a**2 + b**2): # np.sqrt(3^2 + 3^2)
				Z_dist = abs(Zm-Z)
				if Z_dist <= Zmax:
					
					if str(n) not in Mock:
						#print n
						Mock.append(str(n))
						#plt.plot(mock[n][6],mock[n][7],'ro',alpha=0.5)
					Match.append(str(m))
					line+=' match '+str(n+1)+ ' S/N ' + str(mock[n][6]) + ' width ' +str(mock[n][8])
					plt.scatter(float(mock[n][6]),float(mock[n][7]),edgecolor='r',facecolor='none',linewidth=5)
					sn_a.append(mock[n][6])
					w_a.append(mock[n][7])
			
		line+='\n'	

		#if data[m][31] > 3.:
		g.write(line)	
	#plt.legend()

	g.close()
	matches = np.size(Mock)
	print "detected", M,"matched", matches, "multiple matches", np.size(Match)-matches
	#print matches_array
plt.xlabel("S/N")
plt.ylabel("width [km/s]")

plt.show()
