from astropy import units as u
import astropy.coordinates as coord
import numpy as np
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import astropy.coordinates as coord
name,RA,Dec = np.loadtxt("ALMACAL_deep_fields.txt",dtype=str, unpack=True,usecols=(0,1,2))

red='#b70000'
yellow='#ffc600'
blue='#000f9e'
green='#0f8400'
blue2='#31698a'
violet='#551a8b'
data = SkyCoord(RA,Dec,unit=(u.hourangle,u.deg))
plt.rcParams['axes.axisbelow'] = True

ra=coord.Angle(RA,unit=u.hour)
ra_wrap = ra.wrap_at(180*u.degree)
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection="mollweide")

ax.scatter(ra_wrap.radian, data.dec.radian, c=blue)
ax.grid(True)
plt.title("ALMACAL deep")
plt.tight_layout()
fig.savefig('ALMACAL_deep_sky.pdf')
plt.show()
