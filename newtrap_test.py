# test program for newtrap
# Paul H Alfille

import newtrap as nr
import numpy as np

import matplotlib.pyplot as plt
import random

def f( x ):
    return x ** 2 + random.random()

lo = 0
hi = 10
err = .01
target = 4
nr = nr.NewtRap( target, err, lo, hi )

fig = plt.figure()

xs = np.zeros( 200 )
ys = np.zeros( 200 )

y = 0
for i in range(100):
    x = nr.next(y)
    xs[i] = x
    y = f(x)
    ys[i] = y

nr.target = 6    
print(nr.target)
for i in range(100,200):
    x = nr.next(y)
    xs[i] = x
    y = f(x)
    ys[i] = y


plt.plot(xs)
plt.plot(np.clip(ys,lo,hi))
plt.show()

