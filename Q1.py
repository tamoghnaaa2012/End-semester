import math

def N_raphson(f, derivative, x0, tol):
    step = 0
    last_X = x0
    next_X = last_X + 1 * tol  # "different than lastX so loop starts OK
    #print("{:<12} {:<12}".format("Iterations", "Absolute error"))
    while (abs(next_X - last_X) >= tol) and step < 200:  # this is how you terminate the loop - note use of abs()
        new_Y = f(next_X)
        step = step + 1
        last_X = next_X
        next_X = last_X - new_Y / derivative(f, last_X, tol)  # update estimate using N-R
        error = abs(next_X - last_X)

        #file.write("{:^12} {:<12.3e}\n".format(step, error))
        #print("{:^12} {:<12.3e}".format(step, error))

    #file.close()
    return next_X




def f(x):
    return ((x-5)*math.exp(x) +5)


def derivative(f,x,h):
    s = (f(x+h)-f(x-h))/(2*h)
    return s

h=0.001
x=7
ans = N_raphson(f,derivative,x,h)

print ("Root is",ans)     #Ans is the root
#ans = hc/(\lambda)*kT

b = 6.626*pow(10,-34)*3*pow(10,8)/(ans*1.381*pow(10,-23))
print("b = the value of \lambda  is  ",b)



#"C:\Users\TAMOGHNA PATTANAYAK\AppData\Local\Programs\Python\Python39\python.exe" A:/Code/Q1.py
#Root is 4.965114236039297
#b = the value of \lambda  is   0.0028990103282305335

#Process finished with exit code 0
