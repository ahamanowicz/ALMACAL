from astropy import units as u
import astropy.coordinates as coord
import numpy as np
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import matplotlib.patches as patches
import matplotlib
#params = {'text.usetex': False, 'mathtext.fontset': 'stix'}
#plt.rcParams.update(params)

matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['text.latex.unicode'] = True
red='#b70000'
yellow='#ffc600'
blue='#000f9e'
green='#0f8400'
blue2='#31698a'
violet='#551a8b'
grey='#cccccc'
grey_text="#666666"
def redshifted_freq(freq0,z):
	freq = freq0/(1+z)
	return freq

CO=[115.27,230.538,345.796,461.041,576.268,691.473,806.652] #GHz
CO_label = ["CO(1-0)","CO(2-1)","CO(3-2)","CO(4-3)","CO(5-4)","CO(6-5)","CO(7-6)"]
CII = 1900.537 #GHz

Z = np.linspace(0,8,10000)

band3=np.linspace(83,116,1000)#GHz 116-83
b3w = 116-83
band6=np.linspace(211,275,1000)#GHz
b6w = 275-211
band7=np.linspace(275,373,1000)#GHz
b7w=373-275
band8=np.linspace(385,500,1000)#GHz
b8w=500-385
band9=np.linspace(602,720,1000)#GHz
b9w = 720-602

band4 = np.linspace(125,163,1000)#GHz
b4w = 163-125
band5 = np.linspace(163,211,1000)#GHz
b5w = 211-163
n = 8
colors = plt.cm.brg(np.linspace(0,1,n))
angle = 30

fig = plt.figure(figsize=(8,7))
ax = fig.add_subplot(111)
for co,i in zip(CO,range(n-1)):
	
	plt.plot(Z,redshifted_freq(co,Z),color=colors[i], ls='--')


plt.text(0.25,redshifted_freq(CO[0],0.25)+10,CO_label[0], rotation=-15,color=colors[0],fontsize=12,rotation_mode='anchor')
plt.text(0.25,redshifted_freq(CO[1],0.25)+10,CO_label[1], rotation=-30,color=colors[1],fontsize=12,rotation_mode='anchor')
plt.text(0.25,redshifted_freq(CO[2],0.25)+10,CO_label[2], rotation=-40,color=colors[2],fontsize=12,rotation_mode='anchor')
plt.text(0.25,redshifted_freq(CO[3],0.25)+10,CO_label[3], rotation=-50,color=colors[3],fontsize=12,rotation_mode='anchor')
plt.text(0.25,redshifted_freq(CO[4],0.25)+10,CO_label[4], rotation=-55,color=colors[4],fontsize=12,rotation_mode='anchor')
plt.text(0.25,redshifted_freq(CO[5],0.25)+10,CO_label[5], rotation=-60,color=colors[5],fontsize=12,rotation_mode='anchor')
plt.text(0.25,redshifted_freq(CO[6],0.25)+10,CO_label[6], rotation=-60,color=colors[6],fontsize=12,rotation_mode='anchor')



plt.plot(Z,redshifted_freq(CII,Z),color=colors[n-1],ls='--')
plt.text(2.5,redshifted_freq(CII,2.5)+10,"[CII]", rotation=-30,color=colors[n-1],fontsize=12)

ax.add_patch(patches.Rectangle((0,band3[0]),8,b3w,ec=grey, hatch='//////',fc='none'))
ax.text(3.0, np.mean(band3), 'ALMA 3', horizontalalignment='center', verticalalignment='center',fontsize=12,color=grey_text)
ax.add_patch(patches.Rectangle((0,band6[0]),8,b6w,ec=grey, hatch='//////',fc='none'))
ax.text(3.0, np.mean(band6), 'ALMA 6', horizontalalignment='center', verticalalignment='center',fontsize=12,color=grey_text)
ax.add_patch(patches.Rectangle((0,band7[0]),8,b7w,ec=grey, hatch='//////',fc='none'))
ax.text(3.0, np.mean(band7), 'ALMA 7', horizontalalignment='center', verticalalignment='center',fontsize=12,color=grey_text)

ax.add_patch(patches.Rectangle((0,band8[0]),8,b8w,ec=grey, hatch='//////',fc='none'))
ax.text(3.0, np.mean(band8), 'ALMA 8', horizontalalignment='center', verticalalignment='center',fontsize=12,color=grey_text)

ax.add_patch(patches.Rectangle((0,band9[0]),9,b9w,ec=grey, hatch='//////',fc='none'))
ax.text(3.0, np.mean(band9), 'ALMA 9', horizontalalignment='center', verticalalignment='center',fontsize=12,color=grey_text)

"""
ax.add_patch(patches.Rectangle((0,band4[0]),9,b4w,ec=grey, hatch='//////',fc='none'))
ax.text(3.0, np.mean(band4), 'ALMA 4', horizontalalignment='center', verticalalignment='center',fontsize=12,color=grey_text)

ax.add_patch(patches.Rectangle((0,band5[0]),9,b5w,ec=grey, hatch='//////',fc='none'))
ax.text(3.0, np.mean(band5), 'ALMA 5', horizontalalignment='center', verticalalignment='center',fontsize=12,color=grey_text)
"""

plt.xlabel('redshift',fontsize=20)
plt.ylabel(r'$\nu_{obs} $ [GHz] ', fontsize=20)
plt.ylim([0,810])
plt.xlim([0,3.5])
plt.tight_layout()
fig.savefig('CO_ALMA.pdf')
plt.show()
