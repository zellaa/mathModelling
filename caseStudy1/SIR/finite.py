##
import numpy as np
##
# Define initial parameters
finalTime = 10
a = 0
b = 10
finalSpace = b-a
deltaT = 0.01
deltaX = 0.05
J = int(finalSpace/deltaX)
N = int(finalTime/deltaT)
r = 1
a = 1
##
# Placeholder function to generate initial data
def initialiseArrays():
    I = np.zeros((J,N))
    I[30:50,0] = 0.8
    S = 1-I
    R = np.zeros((J,N))
    S = np.reshape(S,(len(S),int(N))).transpose()
    I = np.reshape(I,(len(I),int(N))).transpose()
    R = np.reshape(R,(len(R),int(N))).transpose()
    return S,I,R
##
# Run the explicit Euler scheme
def propagateSIR(S,I,R,N):
    for i in range(0,int(N)-1):
        sUpdate = S[:][i] +((np.roll(S[:][i],1) - 2*S[:][i]+np.roll(S[:][i],-1))*1/(deltaX**2)- a*I[:][i]*R[:][i])*deltaT
        iUpdate = I[:][i] + (r*I[:][i]*S[:][i]- a*I[:][i])*deltaT
        rUpdate = R[:][i] + (a*I[:][i])*deltaT
        S[:][i+1] = sUpdate
        I[:][i+1] = iUpdate
        S[:][i+1] = rUpdate
    return S,I,R
##
S,I,R = initialiseArrays()
S,I,R = propagateSIR(S,I,R,N)