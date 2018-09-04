import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
import aplpy
from matplotlib.patches import Circle, Ellipse
from astropy.wcs import WCS
from astropy import units as u
from mpl_toolkits.axes_grid.anchored_artists import AnchoredAuxTransformBox
import astropy
import sys
########## colours - Hogwart is my Home palette #######
#Hogwart is my Home palette
red='#b70000'
yellow='#ffc600'
blue='#000f9e'
green='#0f8400'

#Fruit salad	
beam='#cccccc'
contour=yellow#'#baffb0'
circ=red#'#ff6363'
########################################
def gaussian(x, amp, cen, wid):
	return amp * exp(-(x-cen)**2 /wid)

calib=sys.argv[1]
cubes= np.loadtxt('fits_list.in', usecols=(0,),dtype='str')
d=20. #defines size of the zoom
with PdfPages(calib+'_candidate_duchamp.pdf') as pdf:
	for cube in cubes:

		freq=np.loadtxt(cube.replace('.fits','_spec.txt'), unpack=True, usecols=(0,))
		freq0,full_w,w_50, w_20,x_peak,y_peak,z_peak,xmin,xmax,ymin,ymax,zmin,zmax=np.loadtxt(cube.replace('.fits','_results.txt'),comments='#',unpack=True, usecols=(9,17,15,16,37,38,39,21,22,23,24,25,26))
		ra,dec = np.loadtxt(cube.replace('.fits','_results.txt'),comments='#',unpack=True, usecols=(5,6), dtype='str')
		detections = np.loadtxt(cube.replace('.fits','.cand'), unpack=True, usecols=(0,),dtype=int)
		nr_det=np.size(freq0)

		if detections != [] :
			for n in detections:
				n=n-1
				print n
				### OPEN files
				obj=np.loadtxt(cube.replace('.fits','_spec.txt'), unpack=True, usecols=(n+1,))  #read the spectrum file for particular cube

				hdul = fits.open(cube) #open cube
				#hdul = fits.open(cube+'.fits') #open cube
				data = hdul[0].data
				wcs = WCS(hdul[0].header)
				wcs=wcs.dropaxis(2)
				wcs=wcs.dropaxis(2)
				if np.size(hdul) > 1.:
					beams=hdul[1].data

					### define beam to plot
					cell= astropy.wcs.utils.proj_plane_pixel_scales(wcs)*3600
					beam_a = np.mean(beams['BMAJ'])/cell[0]
					beam_b = np.mean(beams['BMIN'])/cell[1]
					beam_pa = np.mean(beams['BPA'])
					print "beam = ",beam_a, beam_b, beam_pa
				else: 
					beam_a, beam_b,beam_pa=0.,0.,0.

				#create zer-moment map - collapse the cube over channels with detection
				
				collapsed=data[0,int(zmin[n])]
				
				for z in range(int(zmin[n])+1,int(zmax[n])+1):
					collapsed += data[0,z]

				## spectrum stats
				med=np.median(obj)
				std=np.std(obj)
				mean = np.mean(obj)
				### PLOT

				fig=plt.figure(1,figsize=(15,9))
				#ax=fig.add_subplot(111,projection=wcs) #add sky coords
				ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=2)
				ax2 = plt.subplot2grid((2, 3), (0,2),projection=wcs)
				ax3 = plt.subplot2grid((2,3), (1,0), colspan=3)

				## ax1 - zoomed spectrum over detection
				ax1.step(freq/1e9,obj, '-',c='k', where='mid')
				ax1.axvline(freq0[n]/1e9, linestyle='--', c=blue)
				ax1.axvline((freq0[n]+w_20[n])/1e9,c=red, linestyle='--')
				ax1.axvline((freq0[n]-w_20[n])/1e9,c=red, linestyle='--')
				ax1.axhline(med, c=yellow, linestyle='--')

				ax1.axhline(med+3*std, c=green, linestyle='--')
				ax1.axhline(med-3*std, c=green, linestyle='--')
				ax1.set_xlabel('frequency [GHz]')
				ax1.set_ylabel('flux')
				ax1.set_xlim((freq0[n]-5*w_20[n])/1e9,(freq0[n]+5*w_20[n])/1e9)
				ax1.set_title(cube+" candidate # "+str(n+1))
				#### fits image ######

				# plot the zoomed fits + beam

				box = AnchoredAuxTransformBox(ax2.transData, loc=3,frameon=False,)
				el = Ellipse((0,0), width=beam_b, height=beam_a, angle=beam_pa,fc="none",ec=beam,lw=.7,alpha=0.8) # in data coordinates! hatch='///'
				box.drawing_area.add_artist(el)
				ra = ax2.coords[0]
				ra.set_major_formatter('hh:mm:ss.ss')
				#plt.tick_labels.set_yformat('ddmmss.ss')

				ax2.imshow(collapsed, cmap='gray') 
				sigma_list = np.array([med-7*std,med-6*std,med-5*std,med-4*std,med-3*std,med+3*std,med+4*std,med+5*std,med+6*std,med+7*std])
				ax2.contour(collapsed, colors=contour, levels=sigma_list, alpha=0.5) #levels=np.logspace(-4.7, -3., 10),
				ax2.add_artist(box)
				ax2.set_xlabel('RA')
				ax2.set_ylabel('DEC')
				ax2.set_xlim(x_peak[n]-d,x_peak[n]+d)
				ax2.set_ylim(y_peak[n]-d,y_peak[n]+d)
				#ax2.set_title('integration = '+str(round(time,3))+ ' s')
				lon = ax2.coords[0]
				lat = ax2.coords[1]
				lon.set_ticks(exclude_overlapping=True)
				lat.set_ticks(exclude_overlapping=True)
				lat.set_ticklabel_position('r')
				lat.set_axislabel_position('r')
				#rad = max([abs(x_peak[n]-xmin[n]),abs(x_peak[n]-xmax[n]),abs(y_peak[n]-ymin[n]),abs(y_peak[n]-ymax[n])])
				ax2.plot(x_peak[n],y_peak[n],'+',color=red,ms=7)
				#ax2.add_patch(circ)
				print x_peak[n],y_peak[n]

				#ax2.set_title("RA = "+ra[n]+" DEC = "+dec[n]+"\n S/N = "+str(sn[n]) )

				## ax3 - full spectrum

				ax3.step(freq/1e9,obj, '-',c='k', where='mid')
				ax3.axvline(freq0[n]/1e9, linestyle='--', c=blue)
				ax3.axvline((freq0[n]+w_20[n])/1e9,c=red, linestyle='--')
				ax3.axvline((freq0[n]-w_20[n])/1e9,c=red, linestyle='--')
				ax3.axhline(med, c=yellow, linestyle='--')
				ax3.axhline(med+3*std, c=green, linestyle='--')
				ax3.axhline(med-3*std, c=green, linestyle='--')
				ax3.set_xlabel('frequency [GHz]')
				ax3.set_ylabel('flux')
				#plt.tight_layout()

				#draw the beam
				pdf.savefig()
				plt.close()
		#plt.show()


"""
obj=np.loadtxt(cube+'_spec.txt', unpack=True, usecols=(n+1,))
hdul = fits.open(cube+'.fits')	
hdu = fits.open(cube+'.fits')[0]
wcs = WCS(hdu.header)

data = hdul[0].data
d=15
plt.subplot(projection=wcs)
plt.imshow(data[0,int(z_peak[n])], cmap='gray')
#plt.set_xlabel('pixel')
#plt.set_ylabel('pixel')
#plt.set_xlim(x_peak[n]-d,x_peak[n]+d)
#plt.set_ylim(y_peak[n]-d,y_peak[n]+d)
plt.show()
"""

"""


		obj=np.loadtxt(cube+'_spec.txt', unpack=True, usecols=(n+1,))

		med=np.median(obj)
		std=np.std(obj)
		fig=plt.figure(1,figsize=(15,9))

		ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=2)
		ax2 = plt.subplot2grid((2, 3), (0,2))
		ax3 = plt.subplot2grid((2,3), (1,0), colspan=3)
		ax1.step(freq/1e9,obj, '-',c='k', where='mid')
		ax1.axvline(freq0[n]/1e9, linestyle='--', c=blue)
		ax1.axvline((freq0[n]+w_20[n])/1e9,c=red, linestyle='--')
		ax1.axvline((freq0[n]-w_20[n])/1e9,c=red, linestyle='--')
		ax1.axhline(med, c=yellow, linestyle='--')
		ax1.axhline(med+2*std, c=green, linestyle='--')
		ax1.axhline(med-2*std, c=green, linestyle='--')
		ax1.set_xlabel('frequency [GHz]')
		ax1.set_ylabel('flux')
		ax1.set_xlim((freq0[n]-5*w_20[n])/1e9,(freq0[n]+5*w_20[n])/1e9)

		#### fits image ######

		hdul = fits.open(cube+'.fits')
		
		wcs = WCS(hdul[0].header)

		data = hdul[0].data
		d=15
		ax2.imshow(data[0,int(z_peak[n])], cmap='gray')
		ax2.set_xlabel('pixel')
		ax2.set_ylabel('pixel')
		ax2.set_xlim(x_peak[n]-d,x_peak[n]+d)
		ax2.set_ylim(y_peak[n]-d,y_peak[n]+d)
		rad = max([abs(x_peak[n]-xmin[n]),abs(x_peak[n]-xmax[n]),abs(y_peak[n]-ymin[n]),abs(y_peak[n]-ymax[n])])
		circ = Circle((x_peak[n],y_peak[n]),1.5*rad,ec=red,color='none')
		ax2.add_patch(circ)
		print x_peak[n],y_peak[n]


		ax1.set_title(cube+" candidate # "+str(n+1))
		ax2.set_title("RA = "+ra[n]+" DEC = "+dec[n]+"\n S/N = "+str(sn[n]) )

		ax3.step(freq/1e9,obj, '-',c='k', where='mid')
		ax3.axvline(freq0[n]/1e9, linestyle='--', c=blue)
		ax3.axvline((freq0[n]+w_20[n])/1e9,c=red, linestyle='--')
		ax3.axvline((freq0[n]-w_20[n])/1e9,c=red, linestyle='--')
		ax3.axhline(med, c=yellow, linestyle='--')
		ax3.axhline(med+2*std, c=green, linestyle='--')
		ax3.axhline(med-2*std, c=green, linestyle='--')
		ax3.set_xlabel('frequency [GHz]')
		ax3.set_ylabel('flux')
		plt.tight_layout()

		pdf.savefig()
		plt.close()
"""




