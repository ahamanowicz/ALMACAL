import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as col
import scipy.optimize as so

def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level

data = np.loadtxt("high_mmol.txt")
xdata = data.T[2]

nj = 10 # number of transitions

M = nj * np.size(xdata)
CO_label = ["CO(1-0)","CO(2-1)","CO(3-2)","CO(4-3)","CO(5-4)","CO(6-5)","CO(7-6)", "CO(8-7)", "CO(9-8)", "CO(10-9)"]
CO=[115.27,230.538,345.796,461.041,576.268,691.473,806.652, 921.7997,1036.912393,1151.985452] #GHz

#binning must be the same for all histograms, binedges z: min(z), max(z), and form min(fluxj) to max(fluxj) over all transitrions
nbins=81 #1d number of bins

#find the edges of the ybins 
minj,maxj = [],[]
for j in range(nj):
	if min(np.log10(data.T[j+10])) > -np.inf:
		minj.append(min(np.log10(data.T[j+10])))

	maxj.append(max(np.log10(data.T[j+10])))

minj, maxj = np.array(minj, dtype='float'),np.array(maxj, dtype='float')
print min(xdata), max(xdata),min(minj), max(maxj)
xedges = np.linspace(min(xdata), max(xdata), nbins) ## xbins - range of z
yedges = np.linspace(min(minj), max(maxj), nbins) #xbins - fluxes including all tranistions

#save bins
np.save("xbins",xedges)
np.save("ybins",yedges)


## histogram for each transition
####
Nprob = []
COprob = []
for j in range(nj):
	print j
	z = CO[j]/freq0 -1.
	fig = plt.figure(1,figsize=(8,8))
	ydata = np.log10(data.T[j+10])
	H, xe, ye = np.histogram2d(xdata, ydata, bins=(xedges,yedges))
	print H
	np.save(CO_label[j]+"_hist",H.T)
	

	#f=open("histogram_"+CO_label[j]+".txt", 'w')
	#plt.imshow(H.T, origin='left')
	"""
	if z > 0. and z < 10.:
		print z
		#check which bin you are in for k in range(nbins-1):
		for k in range(nbins-1):
			#print xedges[k],xedges[k+1], yedges[l],yedges[l+1],H.T[k][l]
			x1,x2 = xedges[k],xedges[k+1]
			if x1<=z <=x2: 
				x = k+1
				#print x1,x2,z,k+1
				for l in range(nbins-1):
					y1,y2 = yedges[l],yedges[l+1]
					if y1<=flux <=y2: 
						#print y1,y2,flux,l+1
						y = l+1
						plt.plot(x,y,'or',ms=10)
						plt.plot(y,x,'ob',ms=10)
						Nprob.append(H.T[y][x])
						COprob.append(j)
						print H.T[y][x], j 
				#if y1 <= flu <= y2:
				#	n = hist_array.T[4][k]
				#	prob = n/M		
	#plt.ylim([30,50])
	#plt.xlim([5,25])
	#plt.show()

Nprob = np.array(Nprob, dtype='float')
N = np.sum(Nprob)	
print Nprob
for i in range(np.size(COprob)):
	j = COprob[i]
	z = CO[j]/freq0 -1.
	prob = Nprob[i]/N
	print "transition", CO_label[j], "z =", np.round(z,3), "probability", np.round(prob*100,3)
	#for k in range(nbins-1):
	#	for l in range(nbins-1):
			#print xedges[k],xedges[k+1], yedges[l],yedges[l+1],H.T[k][l]
	#		f.write(str(xedges[l])+", " +str(xedges[l+1])+ ", "+ str(yedges[k]) + ", "+str(yedges[k+1])+", "+str(H.T[k][l])+"\n")
	#f.close()

ydata = np.log10(data.T[10]) #SCO(1-0)
nbins_x,nbins_y=40,40
H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)


x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))
print n, nbins_x, nbins_y
print H, xedges, yedges

fig = plt.figure(1,figsize=(8,8))
H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

pdf = (H*(x_bin_sizes*y_bin_sizes))

low_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.01))
twenty_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.2))
thirty_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.3))
forty_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.4))
fifty_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.5))
one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
eighty_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.8))
ninty_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.9))
two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
levels = [three_sigma, two_sigma, ninty_sigma, eighty_sigma, one_sigma, fifty_sigma, forty_sigma, thirty_sigma, twenty_sigma, low_sigma]

X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
Z = pdf.T


# The viridis colormap is only available since mpl 1.5
extra_args = {}
if tuple(mpl.__version__.split('.')) >= ('1', '5'):
    extra_args['cmap'] = plt.get_cmap('viridis')

plt.contourf(X, Y, Z, levels=levels, origin="lower", alpha=0.75, norm=col.Normalize(vmin=0, vmax=0.01), **extra_args)
plt.colorbar()
plt.xlabel("z", fontsize=15)
#plt.ylabel(r"log10(Mmol) [Msol]", fontsize=15)
#plt.xlim([-1, 11])
# plt.ylim([9.5, 11.5])
fig.savefig("mmol_z.pdf")
plt.show()
"""

