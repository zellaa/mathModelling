##
import numpy as np
##
S = np.array([])
I = np.array([])
R = np.array([])
##
finalTime = 10
a = 0
b = 10
totalSpace = b-a
deltaT = 0.01
deltaX = 0.05
J = int(totalSpace/deltaX)
N = int(finalTime/deltaT)
S = np.array([])
I = np.array([])
R = np.array([])
r = 1
a = 1
##
def initialiseArrays():
    I = np.zeros((J,N))
    I[30:50,0] = 0.8
    S = 1-I
    R = np.zeros((J,N))
    S = np.reshape(S,(len(S),int(N))).transpose()
    I = np.reshape(I,(len(I),int(N))).transpose()
    R = np.reshape(R,(len(R),int(N))).transpose()
    return S,I,R
S,I,R = initialiseArrays()
##
for i in range(0,int(N)-1):
    sUpdate = S[:][i] +((np.roll(S[:][i],1) - 2*S[:][i]+np.roll(S[:][i],-1))*1/(deltaX**2)- a*I[:][i]*R[:][i])*deltaT
    iUpdate = I[:][i] + (r*I[:][i]*S[:][i]- a*I[:][i])*deltaT
    rUpdate = R[:][i] + (a*I[:][i])*deltaT
    S[:][i+1] = sUpdate
    I[:][i+1] = iUpdate
    S[:][i+1] = rUpdate
