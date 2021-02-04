from mylib import *
import math

def func2(x, y, z):
    return (-9.8)

def func1(x, y, z):
    return (z)    #defining d(y)/d(t)=Z

x0 = 0
y0 = 2
zh0 = 40
zl0 = 20

xn = 5
yn = 45
h = 0.05

shoting_method(x0, y0, zh0, zl0, xn, yn, h, func1, func2)
print ("The launch velocity is = the value of Z", )

"""
For boundary condition, Yn =  45
After Langarangian interpolation, The value of Z = 33.25985148514848
And for Z = 33.25985148514848 , The calculated value of Yn = 45.000000000000014
The launch velocity is = the value of Z
"""