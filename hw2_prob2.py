import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d.axes3d import *

def DST(x): #computes discrete since transform of x
    n = np.size(x)
    a = [0 for j in range(0,2*n)]
    a[1:n+1] = x
    a[n+2:2*n+1] = -x[::-1] #odd extension
    y = np.fft.fft(a)   #compute fft
    return np.real(y[1:n+1] / (-2j))

def IDST(x): #computes inverse discrete since transform of x
    n = np.size(x)
    return DST(x) * 2 / (n+1)

def laplace_2d_dirichlet(N): #function for solving -laplacian u = f with Dirichlet BCs.
    h = 1/(N-1)
    x = np.array([j*h for j in range(0,N)])
    y = x
    #evaluate function f.
    F = np.zeros((N,N))
    for i in range(0,N):
        for j in range(0,N):
            #Edit this line for different functions f.
            F[i,j] = np.sin(-10*np.pi*x[i]) * np.sin(5*np.pi*y[j])
            #F[i,j] = x[i]*(x[i]-1)*y[1]*(y[i]-1)
    #compute 2D dst of f.
    Fhat = np.zeros((N,N))
    for i in range(0,N):
        Fhat[i,:] = DST(F[i,:])
    for j in range(0,N):
        Fhat[:,j] = DST(Fhat[:,j])
    #compute coefficients of solution u.
    Uhat = np.zeros((N,N))
    for i in range(0,N):
        for j in range(0,N):
            Uhat[i,j] = Fhat[i,j] / ((i+1)**2 + (j+1)**2)
    Uhat = Uhat / (-np.pi ** 2)
    #compute 2D inverse dst.
    U = np.zeros((N,N))
    for i in range(0,N):
        U[i,:] = IDST(Uhat[i,:])
    for j in range(0,N):
        U[:,j] = IDST(U[:,j])
    return U,F


#Main program

#User input.
N = int(input("Enter value of N."))
U = laplace_2d_dirichlet(N)[0];
F = laplace_2d_dirichlet(N)[1];
#Create plot of u(x).
h = 1/(N-1)
x = np.array([j*h for j in range(0,N)])
#Create X
X = np.tile(x,(N,1))
#Create Y
Y = np.transpose(X)
#plot solution u
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, U, cmap=cm.jet)
plt.show()
#plot f
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, F, cmap=cm.jet)
plt.show()