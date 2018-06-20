import numpy as np
import matplotlib.pyplot as plt
import sys

cubes= np.loadtxt('duchamp_cubes.txt', dtype='str')
f=open('duchamp.in','w')
for cube in cubes:
    
    fig=plt.figure(1,figsize=(6,6))
    stats=imstat(cube,axes=[0,1])['rms']
    chan=np.size(stats)
    med=np.median(stats)
    std=np.std(stats)
    outl=""
    for c in range(chan):
        if abs(stats[c]-med) > std:
            print "outlier = ", c
            if outl=="":
                outl=str(c)
            else:
                outl=outl+','+str(c)
    plt.plot(stats, 'o')
    plt.axhline(med,c='r',linewidth=2)
    plt.axhline(med+std,linestyle='--',c='r',linewidth=2)
    plt.axhline(med-std,linestyle='--',c='r',linewidth=2)
    plt.xlim(-2,chan+2)
    plt.xlabel('channel')
    plt.ylabel('rms')
    plt.title(cube.replace('.image',''))
    plt.tight_layout()
    print outl
    fig.savefig(cube.replace('image','cube')+'_rms.pdf')
    plt.close()
    f.write(cube+' '+outl)
    f.write('\n')
f.close()

