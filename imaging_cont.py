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


uidlist='uid___A001_X1fc_X17.list.J1924+1540'
names = np.loadtxt(uidlist, dtype='str')

###############  parameters for clean ###############

mode='frequency'
interp='linear' #interpolation
weigh='briggs' #weighting
ro=1.0 #robust

################ other parameters #########################

plot_spw = False
min_width = 0.5e9 #minimum cube width [Hz]
max_gap = min_width/2. # max gap between spw merging into cube
field = '0' #there is always only one field in ALMACAL
############################################################

##### Prepare the parameters to create a cube #####

### find the cellsize [in arcsec] ###
for uid in names:
        uid = uid.strip()
        cell=au.pickCellSize(vis=uid,imsize=True,npix=3)[0] #pixelsize in arcsezc
        pix=au.pickCellSize(vis=uid,imsize=True,npix=3)[1][0] # size of an image in pixels

        if cell < 0.5/3. :

                cell_old=cell
                cell = 0.5/3.
                pix_new = int(cell_old/cell *pix)
                pix = au.nextValidImsize(pix_new)

        imagename= uid+"_im"
        clean(vis=uid, imagename=imagename,field=field,mode='mfs',interpolation=interp,imsize=[pix,pix],cell=[cell],weighting=weigh,robust=ro,niter=0,interactive=False,uvtaper=True,outertaper=['0.5arcsec'])

