# test program for newtrap
# Paul H Alfille

import newtrap
import numpy as np

import matplotlib.pyplot as plt
import random

def f( x ):
    return x ** 2 + random.random()

lo = 0
hi = 10
err = .01
target = 4
nr = newtrap.NewtRap( target, err, lo, hi )

fig = plt.figure()

xs = np.zeros( 200 )
ys = np.zeros( 200 )
ts = np.zeros( 200 )

y = 0
for i in range(100):
    x = nr.next(y)
    xs[i] = x
    y = f(x)
    ys[i] = y
    ts[i] = nr.target

nr.target = 6    
for i in range(100,200):
    x = nr.next(y)
    xs[i] = x
    y = f(x)
    ys[i] = y
    ts[i] = nr.target


plt.plot(xs,label="control (x)")
plt.plot(np.clip(ys,lo,hi),label="response (y)")
plt.plot(ts,label="target")
plt.title('NewtRap example -- random error added')
plt.legend()
plt.show()

