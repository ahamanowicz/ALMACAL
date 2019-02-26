############################################################
""" Script for CASApy creating the cubes for ALMACAL
input: list of *.ms files from which you will for a cube

msfile - input file, for clean if the uvcontsub is on would be overwrttten to *.contsub
"""
############################################################

#def fov(freq):
#	c = 3.e8
#	d = 12.
#	FOV = np.degrees(1.22*c/d/freq)*3600 #in arcsec
#	return FOV

##############################################################
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil


uid='uid___A002_Xac3425_X209.ms.split.cal.J1924+1540_B6'
"""
if ".list." in uid: 
                        member_id=uid.split(".")[0]
                        spw = uid.split(".")[4]                
                        calibname_band=uid.split(".")[5]
                        calibname = calibname_band.split("_")[0]
                        band = calibname_band.split("_")[1]
                        imagename= member_id+".cube_"+spw+"."+calibname+'_'+band
else:
                
                        member_id=uid.split(".")[0]
                        spw = uid.split(".")[2]                
                        calibname_band=uid.split(".")[3]
                        calibname = calibname_band.split("_")[0]
                        band = calibname_band.split("_")[1]
                        imagename=uid.replace(".concat.",".cube_")

"""

###############  parameters for clean ###############

mode='frequency'
interp='linear' #interpolation
weigh='briggs' #weighting
ro=1.0 #robust

################ other parameters #########################

plot_spw = False
min_width = 0.5e9 #minimum cube width [Hz]
max_gap = min_width/2. # max gap between spw merging into cube

############################################################

##### Prepare the parameters to create a cube #####

### find the cellsize [in arcsec] ###

cell=au.pickCellSize(vis=uid,imsize=True,npix=3)[0] #pixelsize in arcsezc
pix=au.pickCellSize(vis=uid,imsize=True,npix=3)[1][0] # size of an image in pixels
                
print "before",cell,pix
if cell < 0.5/3. :

        cell_old=cell
        cell = 0.5/3.
        pix_new = int(cell_old/cell *pix)
        pix = au.nextValidImsize(pix_new)


tb.open(uid+'/SPECTRAL_WINDOW')
refFreq = tb.getcol('REF_FREQUENCY')
beam = 1.5*au.primaryBeamArcsec(frequency=refFreq[0]) # beamsize
tb.close()

print "after",cell, pix
print beam, cell*pix
print refFreq


### determine number of fields ###
field = '0' #there is always only one field in ALMACAL


imagename= uid+"_im"

clean(vis=uid, imagename=imagename,field=field,mode='mfs',interpolation=interp,imsize=[pix,pix],cell=[cell],weighting=weigh,robust=ro,niter=0,interactive=False,uvtaper=True,outertaper=['0.5arcsec'])

#rms=imstat(imagename+".image")['rms'][0]

#print imagename
#print 'rms=',rms
#f.write(str(rms)+'\n')


