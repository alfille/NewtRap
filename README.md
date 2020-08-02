# NewtRap
## Newton-Raphson method process control
 The idea is to find the input 'x' that gives an output closest to target

 There is (optional) bounding limits on x: 'lo' and 'hi' 
 you can start with an (optional) initial guess 'x0'
 You can set the allowable error 'error' margin (2-sided)

 Needs no other modules
 (The test program uses matplotlib, random and numpy)

## Method:
 Newton-Raphson for finding a function zero:
 
 x1 = x0 - f(x0)/f'(x0)

 since we don't have a derivative, we need to points to get a difference
 We could choose 2 close points, but that would lead to instability with small output errors
 so we choose the midpoint

## Written by Paul H Alfille 2020

## Usage:
```python
import newtrap
 
 target = 10
 error = .2
 x0 = 1
 lo = 0
 hi = 2
 
 nr = newtrap.NewtRap( target, error=error, lo=lo, hi=hi, x0=x0 )

 x = x0
 
 while True:
    y = my_process(x)
    print(x,y)    
    x = nr.next(y)
```    
## Example from included test program

![Single test plot](test/Figure_1.png)

Here x is chasing a noisy x^2 target

