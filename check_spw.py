import numpy as np
import sys
uid_list = sys.argv[1]

uid_names=np.loadtxt(uid_list, dtype='str',usecols=(0,))

for uid in uid_names:
	
	tb.open(uid+'/SPECTRAL_WINDOW')
	refFreq = tb.getcol('REF_FREQUENCY')
	tb.close()
	tb.open(uid+'/OBSERVATION')
    muid=tb.getcell('PROJECT')
    time = tb.getcell('OBSERVER')
    print "ms", uid, refFreq, time