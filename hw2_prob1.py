import numpy as np
import matplotlib.pyplot as plt

#Function for creating coefficient matrix, A.
#flag = 0 for first problem (not dividing by alpha)
#flag = 1 for first problem (dividing by alpha)
def Abuild(N,alpha,flag):
    h = 1/(N+1) #step size
    a_diag = []
    a_offdiag = []
    for j in range(1,N+1):
        xjplus = (j+0.5)*h
        xjminus = (j-0.5)*h
        ajplus = (np.sin(2*np.pi*xjplus))**2 - alpha * xjplus*(xjplus-1)
        ajminus = (np.sin(2*np.pi*xjminus))**2 - alpha * xjminus*(xjminus-1)
        if (flag==0):
            a_diag.append((ajplus+ajminus))
            a_offdiag.append(-1*ajplus)
        else:
            #modified code for dividing by alpha.
            x = j*h
            ax = (np.sin(2*np.pi*x))**2 - alpha * x*(x-1)
            a_diag.append((ajplus+ajminus)/ax)
            a_offdiag.append(-1*ajplus/ax)

    #remove last element from a_offdiag.
    a_offdiag.pop(-1)
    A = (np.diag(a_offdiag, -1) + np.diag(a_diag, 0) + np.diag(a_offdiag, 1)) / h**2
    return A

#Function for creating vector f(x).
def fbuild(N,alpha,flag):
    h = 1/(N+1) #step size
    f = []
    for j in range(1,N+1):
        x = j*h
        if (flag==0):
            f.append((np.sin(np.pi*x))**2)
        else:
            #modified code for dividing by alpha.
            ax = (np.sin(2*np.pi*x))**2 - alpha * x*(x-1)
            f.append((np.sin(np.pi*x))**2 / ax)
    return f

#Main program

N = int(input("Enter value for N."))
flag = int(input("Enter 0 for first problem or 1 for second problem."))
alpha = float(input("enter value of alpha."))
#Build coefficient matrix, A.
A = Abuild(N,alpha,flag)
#Build RHS vector, f.
f = fbuild(N,alpha,flag)
#Solve linear system for u.
y = np.linalg.solve(A,f)
u = np.zeros(N+2)
u[1:N+1] = y
#Create plot of u(x).
h = 1/(N+1)
x = np.array([j*h for j in range(0,N+2)])
plt.plot(x,u)
plt.legend(['N='+str(N)+', alpha='+str(alpha)])
plt.show()
print(A)
print(f)


