############################################################
""" Script for CASApy creating the cubes for ALMACAL
input: list of *.ms files from which you will for a cube

msfile - input file, for clean if the uvcontsub is on would be overwrttten to *.contsub
"""
############################################################

def fov(freq):
	c = 3.e8
	d = 12.
	FOV = np.degrees(1.22*c/d/freq)*3600 #in arcsec
	return FOV

##############################################################
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil

powers=[]
for n in range(0,15):
        powers.append(2**n)

uid_names=np.loadtxt('names.list', dtype='str',usecols=(0,))

for uid in uid_names:
        
        names=np.loadtxt(uid, dtype='str',usecols=(0,))
        member_id=uid.partition(".list.")[0]
        calibname=uid.partition(".list.")[2]
        print member_id,names[0]

        ###############  parameters for clean ###############

        mode='frequency'
        interp='linear' #interpolation
        weigh='briggs' #weighting
        ro=0.8 #robust

        ################ other parameters #########################

        plot_spw = False
        min_width = 0.5e9 #minimum cube width [Hz]
        max_gap = min_width/2. # max gap between spw merging into cube

        ############################################################

        ##### Prepare the parameters to create a cube #####

        ### find the cellsize [in arcsec] ###
        cell=au.estimateSynthesizedBeam(names[0])/3.
        ### determine number of fields ###
        field = '0' #there is always only one field in ALMACAL

        ##### spectral windows #####
        tb.open(names[0]+'/SPECTRAL_WINDOW')
        refFreq = tb.getcol('REF_FREQUENCY') #read the number of spectral windows
        
        ### determine imsize - take the smallest pixelsize of all spw

        ##### chose spw ####
        spw_boundaries=np.array([])
        spw_index=np.array([]) 	#spw indexes for clean
        #read all spws boundaries, exclude narrow cubes

        for s in range(0,len(refFreq)):
                sp0 = tb.getcell('CHAN_FREQ',rownr=s)[0]
                sp1 = tb.getcell('CHAN_FREQ',rownr=s)[-1]

                if (sp0 > sp1):	sp0,sp1 = sp1,sp0

                width = sp1 - sp0
	
                ### ommit narrow cubes

                if (width >= min_width): 
                        spw_boundaries = np.append(spw_boundaries,[sp0,sp1])
                        spw_index = np.append(spw_index,s)
	
                        if plot_spw == True:
                                #plot resulting spw setup 
                                plt.plot([sp0,sp1],[0,0],'o-')
        del_index = np.array([])

        ### ommit spws repeating frequency coverage	- don't repeat same data

        for i in range(0, len(spw_boundaries)-1,2):
                a,b = spw_boundaries[i],spw_boundaries[i+1]

                for j in range(0, len(spw_boundaries)-1,2):
                        if i!=j:
                                c,d = spw_boundaries[j],spw_boundaries[j+1]
                                if (c >=a and d<=b):
                                        del_index = np.append(del_index, [j,j+1])
			
        spw_boundaries = np.delete(spw_boundaries,del_index)
        spw_index = np.delete(spw_index,del_index[::2]/2.)

        print spw_boundaries
        print spw_index

        ##### creating cubes #####

        ### selecting spws - no big-gap cubes
        spw = []
        taken=np.array([])
        for i in range(0, len(spw_boundaries)-1,2):
                a,b = spw_boundaries[i],spw_boundaries[i+1]
                for j in range(0, len(spw_boundaries)-1,2):
                        if i!=j:
                                c,d = spw_boundaries[j],spw_boundaries[j+1]
                                gap = min(abs(a-d),abs(b-c))
                                if gap < max_gap:
                                        if not spw_index[j/2] in taken:
                                                spw.append(str(int(spw_index[i/2]))+','+str(int(spw_index[j/2])))
                                                taken = np.append(taken,[spw_index[i/2],spw_index[j/2]])
					
        for i in range(0, len(spw_boundaries)-1,2):
                if not spw_index[i/2] in taken:
                        spw.append(str(int(spw_index[i/2])))
                        taken = np.append(taken,[spw_index[j/2]])
        print spw


        uv_spw=[]

        for s in spw:
                if len(s) > 1:
                        for k in s:
                                if k != ',' : uv_spw.append(k)
                else:
                        uv_spw.append(s)
        print uv_spw
        
        ### continuum subvstraction - uvcontsub 
        contsub=[]
        for n in names:
                for s in uv_spw:
                        uvcontsub(vis=n, field='0',fitorder=1,spw=s,fitspw=s)
                        contsub.append(n+'_'+s+'.contsub')
                        shutil.move(n+'.contsub', n+'_'+s+'.contsub')
        
        ### number of pixels for pixsize
        print names[0]
        beam = 1.5*au.primaryBeamArcsec(frequency=refFreq[0])
        print 'b=',beam, 'f=',fov(min(spw_boundaries))
        pix2=int(beam/cell)
        pix = min(filter(lambda x: x >= pix2, powers)) #correct on digestable values
        print "beam=",beam,"pix=",pix


        ### create cubes

        f= open(member_id+'.stat', 'w')
        f.write("beam: "+str(beam)+'\n')
        f.write("pix: "+str(pix)+'\n')
        f.write("n_cubes: "+str(len(spw))+'\n')
        b = names[0].find('.B')
        band=names[0][b-1:]
        f.write("band: "+band+'\n')
        for s in spw:
                if len(s) > 1:
                        
                        fr1=tb.getcell('CHAN_FREQ',rownr=int(s[0]))[0]
                        fr2 = tb.getcell('CHAN_FREQ',rownr=int(s[2]))[-1]
                        
                else:
                        fr1=tb.getcell('CHAN_FREQ',rownr=int(s))[0]
                        fr2=tb.getcell('CHAN_FREQ',rownr=int(s))[-1]
                f.write("frequency_range: "+str(min(fr1,fr2))+", " + str(max(fr1,fr2))+'\n')
                
                images=[]
                for n in contsub:
                        if len(s) > 1:
                                if n[-9] == s[0] or n[-9] == s[2]: images.append(n)
                        else:
                                if n[-9] == s : images.append(n)
                #print images 
                print str(s)
                imagename= member_id+".cube_"+s.replace(',','')+"."+calibname+'_'+band
                clean(vis=images, imagename=imagename,field=field,spw='',mode=mode,interpolation=interp,imsize=[pix,pix],cell=[cell],weighting=weigh,robust=ro,niter=0,interactive=False,uvtaper=True,outertaper=['0.5arcsec'])
                rms=imstat(imagename+'.image')['rms'][0]
                
                print imagename
                print 'rms=',rms
                f.write("rms: "+str(rms)+'\n')
        f.close()

