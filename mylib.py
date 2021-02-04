#  Goss Jordan
#. partial pivot
#  Matrix multiplication
#  Lu decomposition
#  Bisection method
#  First derivative
#  Second derivative
#  Newton raphson method
#  Regula falsi method
#  Laguerre method
#  Synthetic division
#  Mid-point method
#  Trapezoidal method
#  Simpson method
#  Euler method
#  Runge kutta-4
#  Shooting method





########################

def GaussJordan(A, B):
    # the pivot row
    pivot = A[k][k]
    for i in range(k, n):
        A[k][i] /= pivot
    B[k] = B[k] / pivot
    # other rows
    for i in range(n):
        if A[i][k] == 0 or i == k:
            continue
        else:
            term = A[i][k]
            for j in range(k, n):
                A[i][j] = A[i][j] - term * A[k][j]
            B[i] = B[i] - term * B[k]
    return B





def partial_pivot(m, v):
    n = len(m)
    for i in range(n - 1):
        if m[i][i] == 0:
            for j in range(i + 1, n):
                if abs(m[j][i]) > abs(m[i][i]):
                    m[i], m[j] = m[j], m[i]
                    v[i], v[j] = v[j], v[i]

    return (m, v)




def MatrixMultiply(M,A):
    B=[]
    for i in range(len(M)):
        row =[]
        for j in range(len(A[0])):
            row.append(0)
        B.append(row)

    for x in range(len(M)):
        for y in range(len(A[0])):
            for z in range(len(M[0])):
                B[x][y] += M[x][z] * A[z][y]
    return B






def luDecomposition(A, b):
    n = len(A)
    # (1) Extract the b vector
    # b = [0 for i in range(n)]
    # for i in range(0,n):
    #    b[i]=A[i][n]

    L = [[0.0] * n for i in range(n)]
    U = [[0.0] * n for i in range(n)]

    # create the pivot matrix P and the multiplied matrix PA
    partial_pivot(A, b)
    # PA = matrix_mult(P,A)

    # perform the LU decomposition

    for i in range(n):
        for k in range(i, n):
            sum = 0
            for j in range(i):
                sum += (L[i][j] * U[j][k])
            # Evaluating U(i,k)
            U[i][k] = A[i][k] - sum

        for k in range(i, n):
            if (i == k):
                L[i][i] = 1  # diagonal as 1
            else:
                sum = 0
                for j in range(i):
                    sum += (L[k][j] * U[j][i])

                L[k][i] = ((A[k][i] - sum) / U[i][i])

    return (L, U)






def forward_sub(L, b):
    n = len(L)
    y = [[0 * n] for i in range(n)]  # range(start, stop, step)
    for i in range(n):
        sum = 0
        for j in range(i):
            sum = sum + L[i][j] * y[j]

        y[i] = (b[i][0] - sum) / L[i][i]

    return (y)





def back_sub(U, y):
    n = len(U)
    # (6) Perform substitution Ux=y
    x = [0] * n
    for i in range(n - 1, -1, -1):
        for j in range(i + 1, n):
            y[i] = y[i] - U[i][j] * x[j]
        x[i] = y[i] / U[i][i]

    return (x)





def determinant(mat, n):
    temp = [0] * n
    total = 1
    det = 1

    for i in range(0, n):
        index = i

        # finding the index which has non zero value
        while (mat[index][i] == 0 and index < n):
            index += 1

        if (index == n):  # if there is non zero element
            # the determinat of matrix as zero
            continue

        if (index != i):
            # loop for swaping the diagonal element row and index row
            for j in range(0, n):
                mat[index][j], mat[i][j] = mat[i][j], mat[index][j]

                # determinant sign changes when we shift rows
            det = det * int(pow(-1, index - i))

            # storing the values of diagonal row elements
        for j in range(0, n):
            temp[j] = mat[i][j]

            # traversing every row below the diagonal element
        for j in range(i + 1, n):
            num1 = temp[i]  # value of diagonal element
            num2 = mat[j][i]  # value of next row element

            for k in range(0, n):
                mat[j][k] = (num1 * mat[j][k]) - (num2 * temp[k])

            total = total * num1

            # mulitplying the diagonal elements to get determinant
    for i in range(0, n):
        det = det * mat[i][i]

    return int(det / total)  # Det(kA)/k=Det(A);



def firstderv(x,h,func):
    s = (func(x+h)-func(x-h))/(2*h)
    return s
#10
def secndderv(x,h,func):
    s = (func(x+h)+func(x-h)-2*func(x))/(h*h)
    return s


##########
#Assignment 5


def bisection_method(f, a, b, tol):
    step = 0  # iteration counter

    for i in range(1, 200):
        if f(a) * f(b) < 0:
            continue
        else:
            if (abs(f(a)) < abs(f(b))):
                a = a - 0.5 * (b - a)
            else:
                b = b + 0.5 * (b - a)  # upto this portion we are bracketing

    #print("{:<12} {:<12}".format("Iterations", "Absolute error"))  # This is used to tabulate the output
    while (b - a) / 2.0 > tol and step < 200:
        mid_point = (a + b) / 2.0
        # print("Iteration-%d ,mid_point= %f ,f(mid_point) =%f)" %(step,mid_point,f(mid_point)))
        # print('Iteration-%d,')
        if f(mid_point) == 0:
            return (mid_point)  # The midpoint is the root;
        elif f(a) * f(mid_point) < 0:
            b = mid_point

        else:
            a = mid_point

        step = step + 1
        # print(f(a))
        # print(a)
        # print(f(b))
        # print(b)
        # print(f(a)*f(b))
        error = abs(b - a)
        #file.write("{:^12} {:<12.3e}\n".format(step, error))
        #print("{:^12} {:<12.3e}".format(step, error))
    #file.close()
    return (mid_point)




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





def regula_falsi(f, a, b, tol):
    step = 0

    ## a,b is initial guess

    for i in range(1, 200):
        if f(a) * f(b) < 0:
            continue
        else:
            if (abs(f(a)) < abs(f(b))):
                a = a - 0.5 * (b - a)
            else:
                b = b + 0.5 * (b - a)  ## Upto this portion we are bracketing.

    last_c = a
    c = a + 1 * tol
    condition = True
    #print("{:<12} {:<12}".format("Iterations", "Absolute error"))  # This is used to tabulate the output
    while condition and step < 200:
        last_c = c
        c = b - ((b - a) * f(b)) / (f(b) - f(a))
        # print(c)
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
        step = step + 1
        error = abs(c - last_c)
        # print(error)
        condition = abs(f(c)) > tol
        #file.write("{:^12} {:<12.3e}\n".format(step, error))
        #print("{:^12} {:<12.3e}".format(step, error))

    #file.close()
    return (c)


def laguerre(f, f_derivative, s_derivative, alpha, n, tol, file):
    if f(alpha) == 0:
        print("The root is", alpha)
    else:
        alpha1 = alpha + 1 * tol  # "different than lastX so loop starts OK
        step = 0
        print("{:<12} {:<12}".format("Iterations", "Absolute error"))
        while abs(alpha1 - alpha) > tol and step < 200:
            step = step + 1

            G = f_derivative(alpha, 0.4, f) / f(alpha)
            H = G * G - (s_derivative(alpha, 0.4) / f(alpha))
            F = math.sqrt((n - 1) * (n * H - (G * G)))

            if abs(G + F) < abs(G - F):
                a = n / (G - F)
            else:
                a = n / (G + F)
            alpha1 = alpha
            alpha = alpha - a
            error = (alpha1 - alpha)
            file.write("{:^12} {:<12.3e}\n".format(step, error))
            print("{:^12} {:<12.3e}".format(step, error))

            # print(step,error)
            print()
        file.write("{:^12} {:<12}\n".format("Iter number", "error"))
        print("The root is", alpha)

    return (alpha)


def synthetic_div(coeff, alpha):
    while abs(coeff[0]) != 1:
        for k in range(len(coeff)):
            coeff[k] = coeff[k] / coeff[0]

    for l in range(1, len(coeff)):
        coeff[l] = coeff[l] + alpha * coeff[l - 1]
    # print("Coeff",coeff)
    return (coeff)





################
# Assignment 6


def midpoint(f, a, b, n):
    h = float(b - a) / n  # finding midpoint
    result = 0

    for i in range(n):
        result += f((a + h / 2.0))
        a += h
    result *= h  # multiplying with the midpoint
    return result




def trapezoid(f, a, b, n):
    h = (b-a)/float(n) #finding midpoint
    result = f(a) + f(b)

    for i in range(1,n,1):
        result = result + (2*(f(a + i*h))) #finding 2*f(x)
    result *= (h/2.0) #multiplying with dx/2.0

    return result



def simpson(f, a, b, n):
    h = (b-a)/float(n)
    result = f(a) + f(b)

    for i in range(1,n,1):
        if(i%2 ==0):
            result = result + (2*(f(a + i*h)))
        else:
            result = result + (4*(f(a + i*h)))
    result *= (h/3.0)

    return result









################################
#Assignment 7


def euler(x0, y, h, func, x, file):
    temp = 0
    #file.write("{:<15}{:<15}\n".format(x0, y))

    # Iterating till the point at which we
    # need approximation
    while x0 < x:
        temp = y
        y = y + h * func(x0, y)
        x0 = x0 + h

        #file.write("{:<15.6}{:<15.10}\n".format(x0, y))

        # print(  "%.6f"% x0,"      ","%.6f"% y)
    return 0




def rk4(x0, y0, z0, h, func1, func2, x):
    #file.write("{:<15}{:<15}\n".format(x0, y0))
    while x0 < x:
        k1y = h * func1(x0, y0, z0)
        k1z = h * func2(x0, y0, z0)

        k2y = h * func1((x0 + h / 2), (y0 + k1y / 2), (z0 + k1z / 2))
        k2z = h * func2((x0 + h / 2), (y0 + k1y / 2), (z0 + k1z / 2))

        k3y = h * func1((x0 + h / 2), (y0 + k2y / 2), (z0 + k2z / 2))
        k3z = h * func2((x0 + h / 2), (y0 + k2y / 2), (z0 + k2z / 2))

        k4y = h * func1((x0 + h), (y0 + k3y), (z0 + k3z))
        k4z = h * func2((x0 + h), (y0 + k3y), (z0 + k3z))

        y = y0 + (1 / 6 * (k1y + (2 * k2y) + (2 * k3y) + k4y))
        z0 = z0 + (1 / 6 * (k1z + (2 * k2z) + (2 * k3z) + k4z))
        x0 = x0 + h
        y0 = y

        #file.write("{:<15.6}{:<15.10}\n".format(x0, y0))

        # print(  "%.6f"% x0,"      ","%.6f"% y)
        # print( x0, "                ",y)
    return (y0)





def shoting_method(x0, y0, zh0, zl0, xn, yn, h, func1, func2):
    yh = rk4(x0, y0, zh0, h, func1, func2, xn)
    yl = rk4(x0, y0, zl0, h, func1, func2, xn)
    print("For boundary condition, Yn = ", yn)
    if yh > yn and yl < yn:
        while abs(yh - yn) > 0.001 or abs(yl - yn) > 0.001:
            if abs(yh - yn) > abs(yn - yl):
                zh0 = zl0 + (((zh0 - zl0) / (yh - yl)) * (yn - yl))

            elif abs(yh - yn) < abs(yn - yl):
                zl0 = zl0 + (((zh0 - zl0) / (yh - yl)) * (yn - yl))
            yh = rk4(x0, y0, zh0, h, func1, func2, xn)
            yl = rk4(x0, y0, zl0, h, func1, func2, xn)
    elif yh < yn and yl < yn:
        return print("please change your 'zh' ")

    elif yh > yn and yl > yn:
        return print("please change your 'zl' ")

    elif yh < yn and yl > yn:
        return print("please change your 'zl' and 'zh'")
    print("After Langarangian interpolation, The value of Z =", zh0)
    print("And for Z =", zh0, ", The calculated value of Yn =", yh)





def LeastSquare(X,Y):
    n = len(X)
    x1=0;y1=0;x2=0;xy=0;y2=0
    for i in range(n):
        x1+=X[i]
        y1+=Y[i]
        x2+=pow(X[i],2)
        y2 += pow(Y[i], 2)
        xy+=X[i]*Y[i]
    x1=x1/n
    y1=y1/n
    a = ((y1*x2)-(x1*xy))/(x2-(n*pow(x1,2)))
    b = (xy - n*x1*y1)/(x2-n*pow(x1,2))
    Sxx = x2-n*pow(x1,2)
    Syy = y2-n*pow(y1,2)
    Sxy = xy-n*x1*y1
    sigmax = Sxx/n
    sigmay = Syy/n
    covxy = Sxy/n
    r2 = pow(Sxy,2)/(Sxx*Syy)
    return a,b,sigmax,sigmay,covxy,Sxx,Sxy,Syy