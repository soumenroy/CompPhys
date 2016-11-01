from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def spinsys(N, M):
    return 2*np.random.randint(2, size=(N,M)) - 1

def Energy(sys):
    energy = 0.0
    row, colmn = sys.shape
    for i in range(row):
        for j in range(colmn):
            spin = sys[i,j]
            energy += -spin*(sys[(i+1)%row, j] + sys[i,(j+1)%colmn] + sys[(i-1)%row, j] + sys[i,(j-1)%colmn])/2
    return energy




            
def MCSpinArrange(sys, temp):
    row, colmn = sys.shape
    for i in range(row):
        for j in range(colmn):
            spin = sys[i,j]
            Eflip = 2*spin*(sys[(i+1)%N, j] + sys[i,(j+1)%N] + sys[(i-1)%N, j] + sys[i,(j-1)%N])
            rng = np.random.rand()
            boltzman_factor = np.exp(-Eflip/(Kb*temp))
            if Eflip < 0:
                sys[i,j] = -1*spin
            elif  boltzman_factor > rng:
                sys[i,j] = -1*spin
    return sys

Kb = 1.0 # Boltzman Constant    
N, M = 10, 10
num = 30
T = np.linspace(1e-1, 5.0, num)
En = np.zeros(num)
Mt = np.zeros(num)

for i in range(num):
    temp = T[i]
    sys = spinsys(N, M)

    # To equilibrate the system
    steps = 3000
    for j in range(steps):
        sys = MCSpinArrange(sys, temp)
    print i
    # Calculate the physical quantity of the system
    # Physical Qauntity: Energy, Magnetization
    energy = 0.0; mag = 0.0
    psteps = 1000
    for k in range(psteps):
        sys = MCSpinArrange(sys, temp)
        en = Energy(sys)
        energy += en
        Malpha = sys.sum()
        mag += Malpha
        
    En[i] = energy/(psteps*N*M)
    Mt[i] = mag/(psteps*N*M)
    

plt.plot(T, En, 'o', label='Energy')
plt.show()

plt.plot(T, Mt, '*', label='Magnetization')
plt.show()
        
        
    


