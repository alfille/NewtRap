# test program for newtrap
# Paul H Alfille

# version 1m

Study = "Multiple overlay - 3"


import newtrap_12_3 as newtrap

import numpy as np

import matplotlib.pyplot as plt
import random

S = 1
lo = 0*S
hi = 10*S
err = .01*S
target1 = 7*S
target2 = 3*S

def f( x ):
    return 10 - x ** 2 + S*random.random()

#nr = newtrap.NewtRap( target1, err, lo, hi )
nr = newtrap.NewtRap( target1, err, lo, hi )

fig = plt.figure()

xs = np.zeros( 200 )
ys = np.zeros( 200 )
ts = np.zeros( 200 )
ds = np.zeros( 200 )
es = np.zeros( 200 )

first_Pass = True # used for first Label
def Pass():
    global xs,ys,ts,ds,es
    global S,hi,lo,err,target1,target2
    global first_Pass
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
        ds[i] = nr.xhi - nr.xlo # reaching into class to show internals
        es[i] = nr.target - y
        if x == lo:
            print( "LO\t",nr.stage,"X\t",nr.xhi, nr.xmid, nr.xlo, "Y\t",nr.yhi, nr.ymid, nr.ylo ) 
        if x == hi:
            print( "HI\t",nr.stage,"X\t",nr.xhi, nr.xmid, nr.xlo, "Y\t",nr.yhi, nr.ymid, nr.ylo ) 

    plt.subplot(211)
#    plt.plot(xs,".b",label="control (x)" if first_Pass else None) #dot
    plt.plot(xs,"b",label="control (x)" if first_Pass else None) #line
    plt.plot(np.clip(ys,lo,hi),".y",label="response (y)" if first_Pass else None)
    plt.subplot(212)
    plt.plot(np.clip(np.abs(ds),1e-20,1e20),".b",label="dx -- control band" if first_Pass else None)
    plt.plot(np.clip(np.abs(es),1e-20,1e20),".y",label="Absolute response error" if first_Pass else None)

    first_Pass = False

Pass()

plt.subplot(211)
plt.plot(ts,"r",label="target")
plt.title('NewtRap example -- IIR and 3 set')
plt.legend()
plt.subplot(212)
plt.semilogy()

for n in range(100):
    Pass()
    print(len(xs[xs==lo]),'\n')

plt.legend()
plt.show()

