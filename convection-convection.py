#! python3.10

# Temperature distribution on a wall calculated using finite difference method
# The wall is under inflence of convection on both sides

# Import modules

import numpy as np
import matplotlib.pyplot as plt

# Cell count

nx = int(input('Input cell count: '))

# Data

L = 0.04        #Wall thickness
q = 5*10**6     #Heat generation
dx = L/nx       #Infinite part of wall
k = 28          #Conduction coefficient
T0 = 0          #Temp. Left side of wall
Tamb = 30       #Ambient temperature on right side
h = 45          #Convection coefficient

# Matrix creation

A = np.zeros((nx,nx))
B = np.zeros((nx,1))

# Calculation

for i in range(nx):
    if i==0:
        A[i,i]= -((k/dx)+(2*k*h/(h*dx+2*k)))
        A[i,i+1]= k/dx
        B[i,0] = -((2*k*h*Tamb/(h*dx+2*k))+(q*dx))

    if i>0 and i<(nx-1):
        A[i,i]= -2
        A[i,i-1]= 1
        A[i,i+1]= 1
        B[i,0]=-((q*dx**2)/k)

    if i==(nx-1):
        A[i,i]= -((k/dx)+(2*k*h/(h*dx+2*k)))
        A[i,i-1]= k/dx
        B[i,0]=-((2*k*h*Tamb/(h*dx+2*k))+(q*dx))

# Temperature Calculation

T=np.matmul(np.linalg.inv(A),B)
x_num = np.arange((dx/2),(L),dx)

# Plot configs

fig = plt.figure(figsize=(12,8))
plt.plot(x_num,T)
plt.legend([f'Temperature distribution using finite difference method for {nx} cells'])
plt.xlabel('Wall thickness (m)')
plt.ylabel('Temperature (degC)')
plt.show()




