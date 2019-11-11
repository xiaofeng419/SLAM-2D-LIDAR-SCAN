import numpy as np

a = np.arange(9)

a = a.reshape(3,3)
print(a)
c = np.asarray([1, 5, 8])
c = c.reshape(3,1)
b = a[a > c]

a = a.reshape(3,3, 1)
cc = np.asarray(([1, 4])).reshape(1, 1, -1)
dd = cc + a
d = 1