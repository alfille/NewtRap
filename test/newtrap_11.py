#!/bin/python3

# Newton-Raphon method process control
# The idea is to find the input 'x' that gives an output closest to target
#
# There is (optional) bounding limits on x: 'lo' and 'hi' 
# you can start with an (optional) initial guess 'x0'
# You can set the allowable error 'error' margin (2-sided)
#
# needs no other modules
#
# Method:
# Newton-Raphson for finding a function zero:
# x1 = x0 - f(x0)/f'(x0)
#
# since we don't have a derivative, we need to points to get a difference
# We could choose 2 close points, but that would lead to instability with small output errors
# so we choose the midpoint
#
# Written by Paul H Alfille 2020
# MIT license
#
# see https://github.com/alfille/NewtRap
#
# Usage:
# import newtrap
# target = 10
# error = .2
# x0 = 1
# lo = 0
# hi = 2
# nr = newtrap.NewtRap( target, error=error, lo=lo, hi=hi, x0=x0 )
#
# x = x0
# while True:
#     y = my_process(x)
#     print(x,y)
#     x = nr.next(y)

class _slope:
    def __init__(self):
        # keep keep slowly decaying average slope
        self.abs_sum_y = 0
        self.abs_sum_x = 0
        self.sign = -1
    @property
    def slope_average(self):
        if self.abs_sum_x != 0:
            s = self.abs_sum_y / self.abs_sum_x
        else:
            s = 1
        self.sign = -self.sign
        return self.sign * s

    def add_slope(self,dy, dx ):
        if dx != 0:
            self.abs_sum_y = .99999 * self.abs_sum_y + abs(dy)
            self.abs_sum_x = .99999 * self.abs_sum_x + abs(dx)

class _filter():
    def __init__( self, value = None, chain = None ):
        self.value = value
        self.chain = chain
        self._lastx = None
       
    @property
    def lastx(self):
        x = self._lastx
        if self.chain is not None:
            if self.chain.lastx is not None:
                return self.chain.lastx
        return x
    
    @lastx.setter
    def lastx(self,x):
        self._lastx = x
        if self.chain is not None:
            x = self.chain.lastx = x
    
    def get(self, name):
        if name == type(self).__name__:
            return self.value
        elif self.chain:
            return self.chain.get(name)
        else:
            return None

    def set(self, name, value):
        if name == type(self).__name__:
            self.value = value
        if self.chain:
            self.chain.set(name, value )

    def apply( self, x ):
        x = self._apply(x)
        if self.chain:
            x = self.chain.apply( x )
        self._lastx = x
        return x
                
class _lo(_filter):
    # Lower boundary
    def _apply( self, x ):
        if self.value is not None:
            if x < self.value:
                x = self.value
        return x
    
class _hi(_filter):
    # Upper boundary
    def _apply( self, x ):
        if self.value is not None:
            if x > self.value:
                x = self.value
        return x

class _iir(_filter):
    # Infinite filter
    def __init__(self, value=.50, chain = None ):
        # value is decay factor
        super().__init__(value, chain)
        
    def _apply( self, x ):
        if self.value is not None:
            lx = self.lastx
            if lx is not None:
                x = self.value * lx + (1-self.value) * x
        return x

class _end(_filter):
    # does nothing
    def _apply( self, x ):
        return x

class NewtRap():
    # Newton Raphson 's method for control
    # Uses 2 points to find derivative, so needs 2 measurements
    # Remembers internally which measurement
    # Lots of care with all the special cases
    def __init__(self, target=1, error = None, lo=None, hi=None, x0=None):
        self._target = target
        self.slope = _slope()
        
        if error is not None:
            self._error = abs(error)
        elif self._target == 0:
            self._error = .01
        else:
            self._error = .01 * abs(self._target)

        # sort and set lo and hi
        if lo is not None and hi is not None:
            if lo > hi:
                s = hi
                hi = lo
                lo = s

        self.chain = _end()
        if lo is not None:
            self.chain = _lo( value=lo, chain = self.chain )
        if hi is not None:
            self.chain = _hi( value=hi, chain = self.chain )
        self.bound = self.chain # before any filtering
        # IIR filter
        self.chain = _iir( value = .5, chain = self.chain )
        
        # initial x's -- lot's of cases
        if x0 is not None:
            self.xpair0, self.xpair1 = ( x0, 2*x0+1 )
        elif lo is None: # no lo        
            if hi is None:
                self.xpair0, self.xpair1 = ( .5, 1.5 )
            else:
                self.xpair0, self.xpair1 = ( hi, hi-1 )
        elif hi is None: # no lo
            if lo is None:
                self.xpair0, self.xpair1 = ( .5, 1.5 )
            else:
                self.xpair0, self.xpair1 = ( lo, lo+1 )
        else: # bounded
            self.xpair0, self.xpair1 = (lo,hi)
            
        # Set "prior" value
        self.lasty = None
        self.ypair0 = None
        self.ypair1 = None
        
        self.new_settings()
        
    def next( self, value ):

        # prime the pump
        if self.very_first:
            # ignore value (no context)
            self.very_first = False
            self.first = True # which part of the pair?
            self.xpair0 = self.bound.apply(self.xpair0)
            return self.xpair0
            
        # value is from previous x
        y = value - self._target
        
        perfect = (abs(y) <= self._error)
        
        # low pass filter (IIR) y value
        if not perfect:
            if self.lasty is not None:
                alpha = 0
                y = alpha * self.lasty + (1-alpha)*y
            self.lasty = y

        if self.first:
            # from xpair0
            self.ypair0 = y
            if perfect:
                # within tolerances, repeat
                
                # Train filters
                self.chain.apply(self.xpair0)
                self.lasty = y
                
                # (re)send good value
                self.xpair0 = self.bound.apply(self.xpair0)
                return self.xpair0
            else:
                self.first = False
                self.xpair1 = self.chain.apply(self.xpair1)
                return self.xpair1
        else:
            # from xpair1
            self.ypair1 = y
            if perfect:
                # within tolerances, repeat
                
                # Train filters
                self.lasty = y
                self.chain.apply(self.xpair1)
                
                # (re)send good value
                self.xpair1 = self.bound.apply(self.xpair1)
                return self.xpair1
            else:
                self.xpair0, self.xpair1  = self.new_pair()
                self.first = True
                self.xpair0 = self.chain.apply(self.xpair0)
                return self.xpair0

    def new_pair( self ):
        # average and difference
        x1 = .5 * ( self.xpair0 + self.xpair1 )
        y1 = .5 * ( self.ypair0 + self.ypair1 )
        dx = self.xpair0 - self.xpair1
        dy = self.ypair0 - self.ypair1
        
        if dx == 0 or dy == 0:
            # Need an arbitrary angle -- illegal slope otherwise
            if abs(self.ypair0) < abs(self.ypair1):
                x2 = self.xpair0 - self.ypair0/ self.slope.slope_average
            else:
                x2 = self.xpair1 - self.ypair1/ self.slope.slope_average

        else:
            self.slope.add_slope( dy, dx )
            
            # method
            x2 = x1 - y1 * dx / dy
        
        # New bracket
        return ( x2, .5 * (x1 + x2) )

    @property
    def target( self ):
        return self._target
        
    @target.setter
    def target( self, t ):
        # Big jostle
        avg = self.slope.slope_average # avoid second sign reversal
        self.xpair0 += (t - self._target) / avg # deriv is essentially conversion factor
        self.xpair1 += (t - self._target) / avg # deriv is essentially conversion factor
        # New target
        self._target = t
        self.new_settings()

    @property
    def error( self ):
        return self._error
        
    @error.setter
    def error( self, e ):
        self._error = e
        self.new_settings()

    @property
    def lo( self ):
        return self.chain.get("_lo")
        
    @lo.setter
    def lo( self, e ):
        self.chain.set("_lo",e)
        self.new_settings()

    @property
    def hi( self ):
        return self.chain.get("_hi")
        
    @hi.setter
    def hi( self, e ):
        self.chain.set("_hi",e)
        self.new_settings()

    def new_settings( self ):
        # for any change in parameters
        self.very_first = True
        
