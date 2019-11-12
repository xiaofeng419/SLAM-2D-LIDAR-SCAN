import numpy as np
from scipy.ndimage import gaussian_filter


a = np.arange(9)

a = a.reshape(3,3)
print(a)
c = np.asarray([1, 5, 8])
c = c.reshape(3,1)
b = a[a > c]

aa = np.asarray([0, 0,0.0,0,0,1,0,0,0,0,0,0,0])

bb = gaussian_filter(aa, sigma=3)
cc = np.gradient(bb)
bbb = 1