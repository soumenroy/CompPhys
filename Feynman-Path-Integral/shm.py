#!/usr/bin/env python

"""
Ground State Wave Function Of Simple Harmonic Oscillator
"""

import numpy as np
import random, math

#------o--Input Section--o------
n = 33
deltaX = 0.4
epsilon = 0.5      
m = 1000          # Number of Steps

#------o--End Input Section--o----

   
N = np.zeros((n))          # Array for Counter of Point Reselection

#-----o--Energy Function for Particular i th Value--o-----
def e(xi1, xi, xi2, epsilon) :
    e1 = ((2 + epsilon**2)/(2*epsilon))*(xi**2)
    e2 = (1/epsilon)*(xi2 + xi1)*xi
    e = (e1 - e2)/epsilon
    return e

#-----o--New Coordinate Function--o-----
def xkbar(xi1, xi2, epsilon, gamma) :
    mu_k = (xi1 + xi2)/(2 + epsilon**2)
    sigma = math.sqrt(epsilon/(2 + epsilon**2))
    xkbar = mu_k + sigma*gamma
    return xkbar

#---o--Primary Grid Generation--o---
X = np.zeros((n))
for i in range(0,n-1):
        if i <= int(n/4)  :
            X[i] = 0.2 *i
        if i > int(n/4) and i <= int(3*n/4) :
            X[i] = 1.6
        if i > int(3*n/4) and i <= int(n):
            X[i] = 0.2 * (n - i)

for j in range(m):
    E = 0.0; Ebar = 0.0
    alpha = random.uniform(0,1)
    #-----o---Another Primary Grid For Grid Formation--o----
    x = np.zeros((n))
    for i in range(0,n):
        if i <= int(n/4)  :
            x[i] = 0.2 *i
        if i > int(n/4) and i <= int(3*n/4) :
            x[i] = 1.6
        if i > int(3*n/4) and i <= int(n):
            x[i] = 0.2 * (n - i)
    a = int(random.randint(0,int(n-1)))
    gamma = np.random.randn(1)
    for i in range(1,n-1):
        xi = x[i]; xi1 = x[i+1]; xi2 = x[i-1]
        E = E + 0.5*e(x[i+1],x[i],x[i-1],epsilon)
        if i == a : xi = xkbar(xi1, xi2, epsilon, gamma)
             
        if i+1 == a : xi1 = xkbar(xi1, xi2, epsilon, gamma)
            
        if i-2 == a : xi2 = xkbar(xi1, xi2, epsilon, gamma)
            
        Ebar = Ebar + 0.5*e(xi1, xi, xi2, epsilon)
    deltaE = E - Ebar
    #----o--Grid Formation--o-------
    if deltaE < 0 or math.exp(-epsilon*deltaE) > alpha :
        X[a] = xkbar(xi1, xi2, epsilon, gamma)
    if deltaE > 0 and math.exp(-epsilon *deltaE) < alpha :
        X[a] = X[a]
    if X[a] < -3.2:
        X[a] = -3.2
    if X[a] > 3.2:
        X[a] = 3.2
    d = int(((X[a] + 3)*10)/2)
    N[d] = N[d] +1                   # No. Reselection

x0 = -3.2
file = open("shm.dat","w").close()
for l in range(0,n):
    file = open("shm.dat","a")
    file.write(str(x0)+'    '+str(N[l]/m))
    x0 = x0 +0.2
    file.write("\n")
    file.close()

# -----------o--Python Gnuplot--o-----------

import Gnuplot
g = Gnuplot.Gnuplot(debug=1)

g("set term post script")
g("set output 'shm.pdf'")
g("set xlabel 'X_i'")
g("set ylabel '|Psi(x,0)|^2'")
g("set title 'Ground State Wave Function of Harmonic Oscillator' tc lt 8")
g("plot 'shm.dat' using 1:2  w lp lt rgb 'blue' pt 7 notitle")

g("set title 'Fit Plot'")
g("f(x) = a*exp(-b*(x)**2)")
g("fit f(x) 'shm.dat' using 1:2 via a,b")
g("ti = sprintf('a = %.10f  b=%.10f', a, b)")
g("plot [-3.2:3.2][0:0.22] a*exp(-b*(x)**2) ,'shm.dat' using 1:2 w lp lt rgb 'blue' pt 7 t ti")

g("set title 'Histogram Plot'")
g("set boxwidth 0.80 relative")
g("plot 'shm.dat' u 1:2 w boxes lc rgb 'blue' ")
g("set out put")
import os
os.popen("xdg-open shm.pdf")
