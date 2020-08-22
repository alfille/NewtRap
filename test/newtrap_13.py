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

# Add better initial conditions and input filtering

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
        # Needs to alternate
        self.last_IIR_x = None
        
    def _apply( self, x ):
        if self.value is not None:
            lx = self.last_IIR_x
            self.last_IIR_x = x
            if lx is not None:
                x = self.value * x + (1-self.value) * lx
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
        self.chain = _iir( value = 1., chain = self.chain )
        
        # initial x's -- lot's of cases
        if x0 is not None:
            self.xpair0 = self.bound.apply( x0 )
            self.xpair1 = self.bound.apply( .8*x0-1 )
        elif lo is None: # no lo        
            if hi is None: # unbounded
                self.xpair0 = self.bound.apply( 1 )
                self.xpair1 = self.bound.apply( 3 )
            else: # hi only
                self.xpair0 = self.bound.apply( .8*hi - 1 )
                self.xpair1 = self.bound.apply( hi )
        elif hi is None: # lo only
                self.xpair0 = self.bound.apply( lo )
                self.xpair1 = self.bound.apply( 1.2 * lo + 1 )
        else: # bounded
            self.xpair0 = self.bound.apply( .6 * lo + .4 * hi - 1 )
            self.xpair1 = self.bound.apply( .4 * lo + .6 * hi + 1 )
            
        self.inputdecay = .99
        # Set "prior" value
        self.lasty0 = 0
        self.lasty1 = 0
        
        # protect snoopers
        self.ypair0 = None
        self.ypair1 = None
        
        self.new_settings()
        
    def next( self, value ):
        # prime the pump
        if self.very_first:
            # ignore value (no context)
            self.very_first = False
            self.first = True # which part of the pair?
            return self.xpair0
            
        # value is from previous x
        y = value - self._target
        
        perfect = (abs(y) <= self._error)
        
        alpha = self.inputalpha # input filter
        self.chain.set('_iir',alpha)
        if self.first:
            # increasingly weight old values
            self.inputalpha *= self.inputdecay
            # from xpair0
            if perfect:
                # within tolerances, repeat
                self.ypair0 = y
                self.lasty0 = y
                
                # Overtrain filters
                self.chain.apply(self.xpair0)
                
                # (re)send good value
                return self.xpair0
            else:
                # low pass filter (IIR) y value
                self.ypair0 = alpha * y + (1-alpha) * self.lasty0
                self.lasty0 = self.ypair0

                # first half of pair complete, send second half
                self.first = False

                return self.xpair1
        else:
            # from xpair1
            if perfect:
                # within tolerances, repeat
                self.ypair1 = y
                self.lasty = y
                
                # Overtrain filters
                self.chain.apply(self.xpair1)
                
                # (re)send good value
                return self.xpair1
            else:
                # low pass filter (IIR) y value
                self.ypair1 = alpha * y + (1-alpha) * self.lasty1
                self.lasty1 = self.ypair1

                # pair complete with results, calculate next pair
                self.xpair0, self.xpair1  = self.new_pair()
                self.first = True

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
        x2 = self.chain.apply(x2)
        return ( x2, self.bound.apply(.5 * (x1 + x2)) )

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
        self.xpair0 = self.bound.apply( self.xpair0 )
        self.xpair1 = self.bound.apply( self.xpair1 )
        self.new_settings()

    @property
    def hi( self ):
        return self.chain.get("_hi")
        
    @hi.setter
    def hi( self, e ):
        self.chain.set("_hi",e)
        self.xpair0 = self.bound.apply( self.xpair0 )
        self.xpair1 = self.bound.apply( self.xpair1 )
        self.new_settings()

    def new_settings( self ):
        # for any change in parameters
        self.very_first = True
        self.inputalpha = 1.
        
