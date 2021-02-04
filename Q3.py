from mylib import *
import math
import matplotlib.pyplot as plt


def f(a,b,x):
    return b*x+a
def pearson_r(Sxy,Sxx,Syy):
    return ((Sxy)**2 /(Sxx *Syy))


x = [0.00,0.30,0.60,0.90,1.20,1.50,1.80,2.1,2.40,2.70,3.00,3.30]
y = [2.20,1.96,1.72,1.53,1.36,1.22,1.10,1.00,0.86,0.75,0.65,0.60]

a,b,sigmax,sigmay,covxy,Sxx,Sxy,Syy = LeastSquare(x,y)
print("the value of intercept =",a)
print("the value of slope =",b)
print("sigmax =",sigmax)
print("sigmay =",sigmay)
print(covxy)
print("The value of Sxx=",Sxx)
print("The value of Sxy=",Sxy)
print("The value of Syy=",Syy)

ans = pearson_r(Sxy,Sxx,Syy)
print("The value of pearson's r is =",ans)

plt.figure()
plt.plot(x,y)
plt.show()