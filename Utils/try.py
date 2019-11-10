import numpy as np

a = np.arange(9)

a = a.reshape(3,3)
print(a)
c = np.asarray([1, 5, 8])
c = c.reshape(3,1)
b = a[a > c]
print (b)