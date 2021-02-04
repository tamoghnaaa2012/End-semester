from mylib import *


import math
#a is given
a = math.sin(math.pi/8)
#the function
def func(x):
    #The function is given as
    p = 4*pow(9.8,-0.5)*(1/math.sqrt(1-(pow(a,2)*pow(math.sin(x),2))))
    return p
#integrated sum is
sum = simpson(func,0,math.pi/2,10)
print("The time period of oscillation is: ",sum)





#The time period of oscillation is:  2.087320017479592