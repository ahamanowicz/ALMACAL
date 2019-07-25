import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl

# test the code on the bigger sample

#################################################
#Code parameters
mpl.rc ('xtick',labelsize=18)
mpl.rc ('ytick',labelsize=18)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
######
# Take input flux + cenrtral frequency of the b=observed line and calculate the probability for the certain CO transition from Shark SAMS
######p
print "################################################################################ \n Returns probability of the observed line to be a certain CO transition. \n Probabilities are based on the population of galaxies generated from Shark SAMS (Lagos et al 2018)\n for a galaxies with log(Molecular Mass) > 9.5.\n################################################################################"

#freq0 =99.5#float( raw_input("Observed frequency of the line (central in GHz): "))
#flux = np.log10(470.)#np.log10(float(raw_input("Observed integrated flux (in mJy km/s): ")))

#print M
## CO transitions
CO=[115.27,230.538,345.796,461.041,576.268,691.473,806.652, 921.7997,1036.912393,1151.985452] #GHz
CO_label = ["CO(1-0)","CO(2-1)","CO(3-2)","CO(4-3)","CO(5-4)","CO(6-5)","CO(7-6)", "CO(8-7)", "CO(9-8)", "CO(10-9)"]
#J "CO(1-0) - 0 (10) ","CO(2-1) - 1 (11) ","CO(3-2) - 2 (12)","CO(4-3) - 3 (13) ","CO(5-4) -4","CO(6-5) -5 ","CO(7-6) -6", "CO(8-7) -7", "CO(9-8) -8", "CO(10-9) -9"
# in Shark file 10, 11, 12 


#### Create the table of possible transitions 
transition_prob= np.zeros((np.size(CO_label), 4))
#col 0 - J, col1 - z col2 - flux col3 - probability for the transition

##load transition bins
xedges= np.load("xbins.npy")
yedges=np.load("ybins.npy")
nbins=np.size(xedges)

Nprob,COprob,Z=[],[],[]
M=0

f=open("Transition_probability_test.txt", 'w')
f.write("#tested trans, z (Jup, z, Prob)\n")

# read in the sample 
sample = np.loadtxt("test_COtrans.txt", usecols=(0,1,2,3), delimiter=',') #f0, flux [Jykms],z , trans
S = np.size(sample.T[0])
for s in range(S):

	flux = np.log10(sample[s][1]*1e3) #mJy
	freq0 = sample[s][0]
	print "Tested transition z = "+ str(sample[s][2]) + " J = " + CO_label[int(sample[s][3] -1)]
	Nprob,COprob,Z=[],[],[]
	M=0
	f.write(str(sample[s][2])+" "+CO_label[int(sample[s][3] -1)]+" ")
	for j in range(np.size(CO_label)):

		z = CO[j]/freq0 -1.
		flu = np.log10(flux)
		# read in the histogram 
		H = np.load(CO_label[j]+"_hist.npy")
		M += np.sum(H)

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
						if y1<=flux <=y2: 
							#print y1,y2,flux,l+1
							y = l
							#plt.plot(x,y,'or',ms=10)
							#plt.plot(y,x,'ob',ms=10)
							Nprob.append(H[y][x])
							COprob.append(j+1) #Jup
							Z.append(z)
							#print H.T[y][x], j 
		elif z < 0:
			COprob.append(j+1)
			Z.append(str(0))
			Nprob.append(str(0))

	Nprob = np.array(Nprob, dtype='float')
	N = np.sum(Nprob)	
	print Nprob, M
	Mnew=np.sum(Nprob/M)
	Nprob_all = Nprob/M
	for i in range(np.size(COprob)):
		j = COprob[i] - 1 #Jlow
		z = CO[j]/freq0 -1.
		prob = Nprob[i]/M/Mnew
		if z < 0 : z= np.nan
		print "transition", CO_label[j], "z =", np.round(z,3), "probability", np.round(prob*100,3) 
		f.write("(" + str(j+1) +" "+ str(np.round(prob*100,3))+ ")" )
	f.write("\n")

	fig=plt.figure(1, figsize=(7,7))
	plt.plot(COprob,Nprob/M/Mnew, "o--", c='navy',ms=10)
	plt.axvline(sample.T[3][s],c="r")
	plt.xlabel("CO transition Jup",fontsize=18)
	plt.ylabel("Probability",fontsize=18)
	plt.ylim([-0.05, 1.05])
	plt.xlim([-1,11])
	plt.show()




