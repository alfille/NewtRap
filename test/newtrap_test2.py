# test program for newtrap
# Paul H Alfille

# version 2 -- add second plot and clean up loop

import newtrap_1234 as newtrap

import numpy as np

import matplotlib.pyplot as plt
import random

def f( x ):
    return x ** 2 + random.random()

lo = 0
hi = 10
err = .01
target1 = 4
target2 = 6
nr = newtrap.NewtRap( target1, err, lo, hi )

fig = plt.figure()

xs = np.zeros( 200 )
ys = np.zeros( 200 )
ts = np.zeros( 200 )
ds = np.zeros( 200 )
es = np.zeros( 200 )

y = 0
for i in range(200):
    if i == 100:
        nr.target = target2
    x = nr.next(y)
    xs[i] = x
    y = f(x)
    ys[i] = y
    ts[i] = nr.target
    ds[i] = nr.xpair0 - nr.xpair1 # reaching into class to show internals
    es[i] = nr.target - y

plt.subplot(211)
plt.plot(xs,label="control (x)")
plt.plot(np.clip(ys,lo,hi),label="response (y)")
plt.plot(ts,label="target")
plt.title('NewtRap example -- random error added')
plt.legend()
plt.subplot(212)
plt.plot(np.abs(ds),label="dx -- control band")
plt.plot(np.abs(es),label="Absolute response error")
plt.semilogy()
plt.legend()
plt.show()

