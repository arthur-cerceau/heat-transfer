#!python3.10

# Comparison of the temperature distribution inside a wall using finite difference method and analytic method
# The wall is under inflence of convection on one side and isolated on the other
 
# Importing Dependencies
import numpy as np
import matplotlib.pyplot as plt

# Cell Count

nx = int(input('Input cell count: '))

# Data

L = 0.04
q = 5*10**6
dx = L/nx
k = 28
T0 = 0
Tamb = 30
h = 45

# Matrix creation

A = np.zeros((nx,nx))
B = np.zeros((nx,1))

# Calculation

for i in range(0,nx):
    if i==0:
        A[i,i]= -3
        A[i,i+1]= 1
        B[i,0] = -(((q*dx**2)/k)+2*T0)
    if i>0 and i<(nx-1):
        A[i,i]= -2
        A[i,i-1]= 1
        A[i,i+1]= 1
        B[i,0]=-((q*dx**2)/k)
    if i==nx-1:
        A[i,i]= -((k/dx)+(2*k*h/(h*dx+2*k)))
        A[i,i-1]= k/dx
        B[i,0]=-((2*k*h*Tamb/(h*dx+2*k))+(q*dx))
  
# Temperature Calculation

T=np.matmul(np.linalg.inv(A),B)
x_num = np.arange((dx/2),(L),dx)
x = np.arange(0,L,0.001)
c1 = ((q*L+(h*q*L**2/(2*k))-h*T0+h*Tamb))/(h*L+k)
T_ana= -((q*x**2)/(2*k))+c1*x+T0

# Plot configs

fig = plt.figure(figsize=(12,8))
plt.plot(x_num,T)
plt.plot(x,T_ana)
plt.legend([f'Temperature distribution using finite difference method for {nx} cells','Temperature distribution using analytic method'])
plt.xlabel('Wall thickness (m)')
plt.ylabel('Temperature (degC)')
plt.show()