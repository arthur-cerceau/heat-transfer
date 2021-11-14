# Condução 2D em regime estacionário
# Fluxo + Convecção

import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as lin
from pandas import *

# Input data

fluxo = 5000
L = 0.3
W = 0.1
nx,ny = 50,50
dx,dy = L/nx,W/ny
k = 63.9
h = 20
Tinf = 30
A = np.zeros((nx*ny,nx*ny))
B = np.zeros((nx*ny,1))

# Definindo as variáveis a,b,c,d,e,f,g

a = 1/(dx**2)
b = 1/(dy**2)
c = (2*h) / (h*dx**2 + 2*k*dx)
d = (2*h) / (h*dy**2 + 2*k*dy)
e = (2*h*Tinf) / (h*dx**2 + 2*k*dx)
f = (2*h*Tinf) / (h*dy**2 + 2*k*dy)
g = fluxo / (k*dx)

# Calculos
m=-1
for j in range(nx):
    for i in range(ny):
        m+=1
    #Célula 1
        if i==0 and j==0:
            A[m,m] = -(a+b+d)
            A[m,m+1] = a
            A[m,m+nx] = b
            B[m,0] = -(f+g)
    #Células 2
        elif i>0 and i<nx-1 and j==0:
            A[m,m] = -(2*a+b+d)
            A[m,m-1] = a
            A[m,m+1] = a
            A[m,m+nx] = b
            B[m,0] = -f
    #célula 3
        elif i==nx-1 and j==0:
            A[m,m] = -(a+b+c+d)
            A[m,m-1] = a
            A[m,m+nx] = b
            B[m,0] = -(e+f)
    #célula 4
        elif i==0 and j>0 and j<ny-1:
            A[m,m-nx] = b
            A[m,m] = -(a+2*b)
            A[m,m+1] = a
            A[m,m+nx] = b
            B[m,0] = -g
    #célula 5
        elif i>0 and i<nx-1 and j>0 and j<ny-1:
            A[m,m-nx] = b
            A[m,m-1] = a
            A[m,m] = -(2*a+2*b)
            A[m,m+1] = a
            A[m,m+nx] = b
            B[m,0] = 0
    #célula 6
        elif i==nx-1 and j>0 and j<ny-1:
            A[m,m-nx] = b
            A[m,m-1] = a
            A[m,m] = -(a+2*b+c)
            A[m,m+nx] = b
            B[m,0] = -e
    #celula 7
        elif i==0 and j==ny-1:
            A[m,m-nx] = b
            A[m,m] = -(a+b+d)
            A[m,m+1] = a
            B[m,0] = -(f+g)
    #célula 8
        elif i>0 and i<nx-1 and j==ny-1:
            A[m,m-nx] = b
            A[m,m-1] = a
            A[m,m] = -(2*a+b+d)
            A[m,m+1] = a
            B[m,0] = -f
    #celula 9
        elif i==nx-1 and j==ny-1:
            A[m,m-nx] = b
            A[m,m-1] = a
            A[m,m] = -(a+b+c+d)
            B[m,0] = -(e+f)
            

T=np.linalg.solve(A,B)
TT=np.zeros((nx,ny))
mm=-1
for ii in range(ny):
    for jj in range(nx):
        mm+=1
        TT[ii,jj]=T[mm]


fig = plt.figure(figsize=(12,8))
x = np.arange(dx/2,L,dx)
y = np.arange(dy/2,W,dy)
plt.contourf(x,y,TT,cmap='gist_rainbow_r',levels=1000)
plt.xlabel('Wall lenght (m)')
plt.ylabel('Wall width (m)')
plt.colorbar(label='Temperature (degC)')
plt.show()
