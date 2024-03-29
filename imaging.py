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

fieldsfile = 'fields_extra.list'

g = open("/scigarfs/opsw/work/ahamanow/ALMACAL/deep/"+fieldsfile+".log",'w')
calibname = np.loadtxt(fieldsfile,dtype='str') #run over all field names
QsoList=np.loadtxt("qsoclean.list",dtype='str')
for c in calibname:
        c = c.split('.')[0]
        #enter the calibrator folder
        g.write("### calibrator "+ c + "\n")
        os.chdir(c) # enter the calibrator folder
	print "field", c
        #files with the imaging parameters for all cubes in the calibrator folder
        f = open(c+'.stat', 'w') 
        f.write("#cubename\t beam\t pix\t cell\t \n")

        #create the concat.list
        concat_list = os.listdir('.')
        uid_names = [b for b in concat_list if ".concat" in b]
        N = len(uid_names)
        print N,"cubes"
        #open the file with conac ms files names
        #uid_names=np.loadtxt('concat.list', dtype='str',usecols=(0,))
        #N = len(uid_names)
	#print uid_names
        for m in range(N):
	       
                uid=uid_names[m]
                print str(m) + " \ " + str(N)
                g.write(uid+'\n')
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
                concatname=str(uid)
                print member_id,calibname, concatname
                #print spw,band, uid
                ###############  parameters for clean ###############

                mode='frequency'
                interp='linear' #interpolation
                weigh='briggs' #weighting
                ro=1.0 #robust
                niter=0
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
                ######  stats file  ######
                
                f.write(imagename+".image"+"\t"+ str(beam)+"\t"+str(pix)+'\t'+str(cell)+'\t')
		f.write("\n")
                ### determine number of fields ###
                field = '0' #there is always only one field in ALMACAL
                #if imagename in QsoList: 
                #        print "QSO res - CLEAN"
                #        niter=100              
                print imagename,concatname
        
                #imagename= member_id+".cube_"+s.rep+"."+calibname+'_'+band+".test_npix5"
                concatname= uid #member_id+".concat."+spw+"."+calibname+'_'+band
                #concat(vis=images,concatvis=concatname,dirtol='20arcsec')
                g.write("#before clean "+imagename+"\n")
                clean(vis=concatname, imagename=imagename,field=field,spw='',mode=mode,interpolation=interp,imsize=[pix,pix],cell=[cell],weighting=weigh,robust=ro,niter=niter,,interactive=False,uvtaper=True,outertaper=['0.5arcsec'])
                g.write("#after clean "+imagename+"\n")
                #rms=imstat(imagename+".image")['rms'][0]
                        
                #print imagename
                #print 'rms=',rms

                #f.write(str(rms)+'\n')
        
        f.close()
        os.chdir('../')
g.close()
  #from here the is big skip compared to regular cubeformcode, we dont choose spw sice the were chosen already andn we cuse concatenated data

                ##### spectral windows #####
                #tb.open(names[0]+'/SPECTRAL_WINDOW')
                #refFreq = tb.getcol('REF_FREQUENCY') #read the number of spectral windows
                
                ##### chose spw #### [no need for repeating imaging - spw included in the concatenated file ]
                #spw_boundaries=np.array([])
                #spw_index=np.array([])         #spw indexes for clean
                #read all spws boundaries, exclude narrow cubes

                #for s in range(0,len(refFreq)):
                #        sp0 = tb.getcell('CHAN_FREQ',rownr=s)[0]
                #        sp1 = tb.getcell('CHAN_FREQ',rownr=s)[-1]
                #        if (sp0 > sp1):        sp0,sp1 = sp1,sp0
                #        width = sp1 - sp0
                
                #        ### ommit narrow cubes

                #        if (width >= min_width): 
                #                spw_boundaries = np.append(spw_boundaries,[sp0,sp1])
                #                spw_index = np.append(spw_index,s)
                
                #                if plot_spw == True:
                                        #plot resulting spw setup 
                #                        plt.plot([sp0,sp1],[0,0],'o-')
                #del_index = np.array([])

                ### ommit spws repeating frequency coverage - don't repeat same data

                #for i in range(0, len(spw_boundaries)-1,2):
                #        a,b = spw_boundaries[i],spw_boundaries[i+1]

                #        for j in range(0, len(spw_boundaries)-1,2):
                #                if i!=j:
                #                        c,d = spw_boundaries[j],spw_boundaries[j+1]
                #                        if (c >=a and d<=b):
                #                                del_index = np.append(del_index, [j,j+1])
                                
                #spw_boundaries = np.delete(spw_boundaries,del_index)
                #spw_index = np.delete(spw_index,del_index[::2]/2.)

                #print spw_boundaries
                #print spw_index

                ##### creating cubes #####

                ### selecting spws - no big-gap cubes [ no need, spw already chosen in the concatenated file]
                #spw = []
                #taken=np.array([])
                #for i in range(0, len(spw_boundaries)-1,2):
                #        a,b = spw_boundaries[i],spw_boundaries[i+1]
                #        for j in range(0, len(spw_boundaries)-1,2):
                #                if i!=j:
                #                        c,d = spw_boundaries[j],spw_boundaries[j+1]
                #                        gap = min(abs(a-d),abs(b-c))
                #                        if gap < max_gap:
                #                                if not spw_index[j/2] in taken:
                #                                        spw.append(str(int(spw_index[i/2]))+','+str(int(spw_index[j/2])))
                #                                        taken = np.append(taken,[spw_index[i/2],spw_index[j/2]])
                                                
                #for i in range(0, len(spw_boundaries)-1,2):
                #        if not spw_index[i/2] in taken:
                #                spw.append(str(int(spw_index[i/2])))
                #                taken = np.append(taken,[spw_index[j/2]])
                #print spw
                #uv_spw=[] # all files go first through uvcontsub only then are concatenated - skip when start from concatenated

                #for s in spw:
                #        if len(s) > 1:
                #                for k in s:
                #                        if k != ',' : uv_spw.append(k)
                #        else:
                #                uv_spw.append(s)
                #print uv_spw
                
                ### continuum subbstraction - uvcontsub  [skip when using the concatenated cubes]
                
                #contsub=[]
                #for n in names:
                #        for s in uv_spw:
                #                uvcontsub(vis=n, field='0',fitorder=1,spw=s,fitspw=s)
                #                contsub.append(n+'_'+s+'.contsub')
                #               shutil.move(n+'.contsub', n+'_'+s+'.contsub')
                
                ######  stats file  ######
                #b = names[0].find('.B')
                #band=names[0][b-1:]
                
                #f.write(member_id+".cube_"+s.replace(',','')+"."+calibname+'_'+band+".image\t"+ str(beam)+"\t"+str(pix)+'\t'+str(cell)+'\n')
                #tb.open(names[0]+'/SPECTRAL_WINDOW')
                #for s in spw:
                #        if len(s) > 1:
                #                
                #                fr1=tb.getcell('CHAN_FREQ',rownr=int(s[0]))[0]
                #                fr2 = tb.getcell('CHAN_FREQ',rownr=int(s[2]))[-1]
                #                
                #        else:
                #                fr1=tb.getcell('CHAN_FREQ',rownr=int(s))[0]
                #                fr2=tb.getcell('CHAN_FREQ',rownr=int(s))[-1]
                
                #        f.write("frequency_range: "+str(min(fr1,fr2))+", " + str(max(fr1,fr2))+'\n')
                        
                #tb.close()

                #for s in spw:
                #        images=[]
                #        for n in contsub:
                #                if len(s) > 1:
                #                        if n[-9] == s[0] or n[-9] == s[2]: images.append(n)
                #                else:
                #                        if n[-9] == s : images.append(n)
                        #print images 
                #        print str(s)

