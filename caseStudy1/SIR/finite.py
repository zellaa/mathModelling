##
import numpy as np
import matplotlib.pyplot as plt
##
# Define initial parameters

beta, gamma = 0.15, 0.04
rho = gamma/beta
deathRate = 30/100000
alpha = deathRate/beta

N0 = 100000
i0, r0 = 100, 0
s0 = N0-i0-r0
s0Hat, i0Hat, r0Hat, d0Hat = (s0/N0), (i0/N0), (r0/N0), 0
y0Hat = s0Hat, i0Hat, r0Hat, d0Hat
##

#finalTime = 10
finalTime = 20
a = 0
b = 2
finalSpace = b-a
deltaT = 0.001
deltaX = 0.01
J = int(round(finalSpace/deltaX))
N = int(round(finalTime/deltaT))
r = 0.02
p = 0.5

##
# Placeholder function to generate initial data
def initialiseArrays():
    S = np.zeros((J,N))
    I = np.zeros((J,N))
    R = np.zeros((J,N))

    S = np.reshape(S,(len(S),int(N))).transpose()
    I = np.reshape(I,(len(I),int(N))).transpose()
    R = np.reshape(R,(len(R),int(N))).transpose()

    S[0] = s0Hat
    I[0] = i0Hat
    R[0] = r0Hat
    return S,I,R
##
# Run the explicit Euler scheme
def propagateSIR(S,I,R,N):
    for i in range(0,int(N)-1):
        sUpdate = S[:][i] +(
                p*(np.roll(S[:][i],1) - 2*S[:][i]+np.roll(S[:][i],-1))*1/(deltaX**2)
                            - I[:][i]*S[:][i])*deltaT
        iUpdate = I[:][i] + (I[:][i]*S[:][i]- rho*I[:][i])*deltaT
        rUpdate = R[:][i] + (rho*I[:][i])*deltaT
        S[:][i+1] = sUpdate
        I[:][i+1] = iUpdate
        R[:][i+1] = rUpdate
    return S,I,R
##
S,I,R = initialiseArrays()
##
S,I,R = propagateSIR(S,I,R,N)
##
tValues = np.linspace(0,finalTime,num=N)/beta
xValues = np.linspace(0,finalSpace,num=J)
X,T = np.meshgrid(xValues,tValues)
fig = plt.figure(figsize=plt.figaspect(0.35))
ax1 = fig.add_subplot(1,3,1,projection='3d')
ax1.contour3D(X,T,S*N0,50)
ax1.set_xlabel('Position in Space')
ax1.set_ylabel('Time')
ax1.set_zlabel('Susceptibles')
ax1.view_init(20, 135)

ax2 = fig.add_subplot(1,3,2,projection='3d')
ax2.contour3D(X,T,I*N0,50)
ax2.view_init(20, 135)
ax2.set_xlabel('Position in Space')
ax2.set_ylabel('Time')
ax2.set_zlabel('Infected')

ax3 = fig.add_subplot(1,3,3,projection='3d')
ax3.contour3D(X,T,R*N0,50)
ax3.view_init(20, 135)
ax3.set_xlabel('Position in Space')
ax3.set_ylabel('Time')
ax3.set_zlabel('Recovered')
plt.show()