# test program for newtrap
# Paul H Alfille

# version 1m

Study = "Multiple overlay"


import newtrap_9 as newtrap

import numpy as np

import matplotlib.pyplot as plt
import random

S = 1
lo = 0*S
hi = 10*S
err = .01*S
target1 = 4*S
target2 = 6*S

def f( x ):
    return x ** 2 + S*random.random()

#nr = newtrap.NewtRap( target1, err, lo, hi )
nr = newtrap.NewtRap( target1, err, lo, hi )

fig = plt.figure()

xs = np.zeros( 200 )
ys = np.zeros( 200 )
ts = np.zeros( 200 )
ds = np.zeros( 200 )
es = np.zeros( 200 )

def Pass():
    global xs,ys,ts,ds,es
    global S,hi,lo,err,target1,target2
    y = 0
    nr.target = target1
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
    plt.plot(xs,".b")
    plt.plot(np.clip(ys,lo,hi),".y")
    plt.subplot(212)
    plt.plot(np.clip(np.abs(ds),1e20,1e-20),".b")
    plt.plot(np.clip(np.abs(es),1e20,1e-20),".y")

Pass()

plt.subplot(211)
plt.plot(xs,".b",label="control (x)")
plt.plot(np.clip(ys,lo,hi),".y",label="response (y)")
plt.plot(ts,"r",label="target")
plt.title('NewtRap example -- filter input')
plt.legend()
plt.subplot(212)
plt.plot(np.abs(ds),".b",label="dx -- control band")
plt.plot(np.abs(es),".y",label="Absolute response error")
plt.semilogy()

for n in range(100):
    Pass()

plt.legend()
plt.show()

