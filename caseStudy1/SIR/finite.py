##
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
##
# Define initial parameters
beta, gamma = 0.15, 0.04
rho = gamma/beta
deathRate = 30/100000
alpha = deathRate/beta
# N0 is the initial population
N0 = 100000
i0, r0 = 100, 0
s0 = N0-i0-r0
s0Hat, i0Hat, r0Hat, d0Hat = (s0/N0), (i0/N0), (r0/N0), 0
y0Hat = s0Hat, i0Hat, r0Hat, d0Hat
finalTime = 80
x0 = 0
xJ = 4
finalSpace = xJ-x0
deltaT = 0.001
deltaX = 0.01
# Final lengths of arrays
J = int(round(finalSpace/deltaX))
N = int(round(finalTime/deltaT))
diffusionScaling = 0.05
##
# Placeholder function to generate initial data
def initialiseArrays():
    S = np.zeros((J,N))
    I = np.zeros((J,N))
    R = np.zeros((J,N))
    D = np.zeros((J,N))

    S = np.reshape(S,(len(S),int(N))).transpose()
    I = np.reshape(I,(len(I),int(N))).transpose()
    R = np.reshape(R,(len(R),int(N))).transpose()
    D = np.reshape(D,(len(D),int(N))).transpose()

    S[0][100:200] = s0Hat
    I[0] = i0Hat
    R[0] = r0Hat
    D[0] = d0Hat
    return S,I,R,D
##
# Run the explicit Euler scheme
def propagateSIR(S,I,R,D,N):
    for i in range(0,int(N)-1):
        sUpdate = S[:][i] +(
                diffusionScaling*(np.roll(S[:][i],1) - 2*S[:][i]+np.roll(S[:][i],-1))*1/(deltaX**2)
                            - I[:][i]*S[:][i])*deltaT
        iUpdate = I[:][i] + (I[:][i]*S[:][i]- rho*I[:][i]-alpha*I[:][i])*deltaT
        rUpdate = R[:][i] + (rho*I[:][i])*deltaT
        dUpdate = D[:][i] + (alpha*I[:][i])*deltaT
        S[:][i+1] = sUpdate
        I[:][i+1] = iUpdate
        R[:][i+1] = rUpdate
        D[:][i+1] = dUpdate
    return S,I,R,D
##
S,I,R,D = initialiseArrays()
##
S,I,R,D = propagateSIR(S,I,R,D,N)
##
def sirPlotter(S,I,R,D,N,J,finalTime,finalSpace,beta):
    tValues = np.linspace(0,finalTime,num=N)/beta
    xValues = np.linspace(0,finalSpace,num=J)
    X,T = np.meshgrid(xValues,tValues)
    fig = plt.figure(figsize=(5,5))
    ax1 = fig.add_subplot(2,2,1,projection='3d')
    ax1.plot_surface(X,T,S,cmap=cm.coolwarm)
    ax1.set_xlabel('Position in Space')
    ax1.set_ylabel('Time (days)')
    ax1.set_zlabel('Susceptibles')
    ax1.view_init(20, 135)

    ax2 = fig.add_subplot(2,2,2,projection='3d')
    ax2.plot_surface(X,T,I,cmap=cm.coolwarm)
    ax2.view_init(20, 135)
    ax2.set_xlabel('Position in Space')
    ax2.set_ylabel('Time (days)')
    ax2.set_zlabel('Infected')

    ax3 = fig.add_subplot(2,2,3,projection='3d')
    ax3.plot_surface(X,T,R,cmap=cm.coolwarm)
    ax3.view_init(20, 135)
    ax3.set_xlabel('Position in Space')
    ax3.set_ylabel('Time (days)')
    ax3.set_zlabel('Recovered')

    ax4 = fig.add_subplot(2,2,4,projection='3d')
    ax4.plot_surface(X,T,D,cmap=cm.coolwarm)
    ax4.view_init(20, 135)
    ax4.set_xlabel('Position in Space')
    ax4.set_ylabel('Time (days)')
    ax4.set_zlabel('Dead')
    plt.suptitle('Plot of S,I,R with a diffusivity model (percentage of populace)')
    plt.show()

sirPlotter(S,I,R,D,N,J,finalTime,finalSpace,beta)
