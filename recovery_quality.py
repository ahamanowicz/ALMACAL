# how well we retireve the mock detections parametersi
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rc ('xtick',labelsize=15)
mpl.rc ('ytick',labelsize=15)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

red='#d62d20'
yellow='#ffa700'
blue='#0057e7'
green='#008744'
orange='#f37735'

cubes = np.loadtxt("fits_list_all.in", usecols=(0,), dtype='str')
SN_frac, W_frac, SN_in, W_in, F_frac, F_in=[],[],[],[],[],[]
for c in cubes:
	c2 = c.replace(".fits", ".matches")
	### SN width from the matches - S/N and width in channel retrieved
	sn_in, w_in, sn_out, w_out, f_in, f_out = np.loadtxt(c2, usecols=(16,18,8,10,20,12), unpack=True, dtype='float')
	SN_frac = np.append(SN_frac, sn_out/sn_in)
	W_frac = np.append(W_frac, w_out/w_in)
	SN_in = np.append(SN_in, sn_in)
	W_in = np.append(W_in, w_in)
	F_frac = np.append(F_frac, f_out/f_in)
	F_in = np.append(F_in, f_in)

fig=plt.figure(1,figsize=(12,6))
plt.subplot(121)
plt.scatter(SN_in,SN_frac, c=red )
plt.ylabel(r"S/N$_{in}$ / S/N$_{out}$", fontsize=18)
plt.xlabel(r"S/N$_{in}$", fontsize=18)
plt.axhline(1, c='k')
#plt.subplot(132)
#plt.scatter(F_in,F_frac, c=green)
#plt.ylabel(r"Peak FLux$_{in}$ / Peak Flux$_{out}$", fontsize=18)
#plt.xlabel(r"Peak Flux$_{in}$ [mJy/beam]", fontsize=18)
#plt.xlim([0])
#plt.axhline(1, c='k')
plt.subplot(122)
plt.scatter(W_in,W_frac, c=blue )
plt.ylabel(r"Width$_{in}$ / Width$_{out}$", fontsize=18)
plt.xlabel(r"Width$_{in}$ [channels]", fontsize=18)
plt.axhline(1, c='k')
fig.savefig("recovery.pdf")
plt.show()