import numpy as np
from scipy.ndimage import gaussian_filter


a = np.arange(9)

a = a.reshape(3,3)
print(a)

xx = np.array(([0,1], [1,2]))
yy = np.array(([1,2], [1,1]))

a[xx, yy] = 10
print (a)