import matplotlib.pyplot as plt
import numpy as np

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

#for mu, sig in [(-1, 1), (0, 2), (2, 3)]:
#    plt.plot(gaussian(np.linspace(-3, 3, 120), mu, sig))	
plt.figure(1)
plt.plot(gaussian(np.linspace(-3, 3, 120), 1,1))
plt.plot()