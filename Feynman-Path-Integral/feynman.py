import numpy as np
import math
import cmath
import trace
from numpy import linalg as la
m = 1.0
hbar = 1.0
xmax =2.0 #input("entr the value of 'xmax' : ")
xmin =-2.0 #input("entr the value of 'xmin' : ")
n =200# input("entr the value of 'n' : ")
tmax =10.0# input("enter the value of 'tmax' : ")
tmin =0.0 #input("enter the value of 'tmin' : ")
p =100 #input("enter the value of 'p' : ")
x = np.linspace(xmax,xmin,n)
t = np.linspace(tmax,tmin,p)
deltaT = ((tmax-tmin)/p)
deltaX = (xmax-xmin)/n
delta = deltaX**2
DT = 512
# This parameter determines the expansion of the procedure
Nt = 3*p
# The constant C appear from propagator

pp = (2*math.pi*hbar*deltaT)
C = cmath.sqrt(-1j/pp)
print C

#Pre -allocation of matrix (short time propagator)

K = np.zeros((n,n),dtype=complex)

#Pre -allocation of matrix (short time propagator)

U = np.zeros(((p+Nt),1),dtype=complex)

#Pre -allocation of vector (for plotting purposes)

r = np.zeros((n,1))

#Pre -allocation of vector (for plotting purposes)

s = np.zeros((n,1))

#Pre -allocation of vector (for discrete fourier transformation)

dft = np.zeros(((DT*p),1))



for j in range(n):
    for k in range(n):
        L = 0.5*((x[j]-x[k])/deltaT)**2 + 0.5*((x[j]+x[k])/2)**2
        #print L
        CC =L*C*cmath.exp(1j*deltaT/hbar)
        #print CC
        K[j,k] = CC
        #print K[j,k]
G = np.array(K, dtype=complex)
print G
# parameter tmp to reduce number of matrix multiplication

tmp = K*deltaX

#trace of the first finite time propagator (G)
Tr = 0
for s in range(n):
    Tr = Tr + G[s,s]
    print G[s,s]
U[1]=Tr
print U[1]
trac = np.array(G)
trac = np.trace(G)
print trac

"""
MATRIX CLASS:
mymatrix = [[1,2],[3,4]]
from numpy import *
from numpy.linalg import *
a = mymatrix
trace = trace(a)
inverse = inverse(a)
a = array(a) # to show like matrix
print a
eigenvalues = eig(a)
unit = eye(3) # to produce 3*3 unit matrix
transpose = a.transpose()
matrix_product = dot(a,a)
b = array([1,2])
x = solve(a,b) # solution of symstem of equation
from numpy import linalg as la
norm = la.norm(a)
"""

for recur in range (2,p):
    
    #calculation of finite time propagator as the product of short time propagator

    G = G*tmp
    # normalizing the propagator (thus U will be unitary)
    
    norm_G = la.norm(G*deltaX)
    G = G/norm_G

    # one eighth through the time interval
    
    if recur==round(p/16):
        for m in range (1,n):
            SUM = 0
            for  j in range (1,n):
                for  k in range (1,n):

                    # evaluating the normalization at given time

                    SUM = SUM + np.conjugate(G[j,k])*G[j,m]*delta
                    print SUM
        r[m] = SUM.real
        s[m] = SUM.imag

