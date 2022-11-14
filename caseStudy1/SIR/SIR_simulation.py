import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

default = {'beta':0.15,
            'gamma':0.04,
            'initialState':[1e6-100,100,0],
            'deathRate':3e-4,
            'T':365,
            'dt':1}

class SIR:
    def __init__(self,parameters=default):
        print('initialising')
        self.beta = parameters['beta']
        self.gamma = parameters['gamma']
        self.deathRate = parameters['deathRate']
        self.rho = self.gamma/self.beta
        self.alpha = self.deathRate/self.beta
        self.T = parameters['T']
        self.dt = parameters['dt']
        self.ts = np.arange(0,self.T,self.dt)
        self.s0,self.i0,self.r0 = parameters['initialState']
        self.N = np.sum(parameters['initialState'])
        self.nonDimenionalise()
    
    def nonDimenionalise(self):
        self.s0Hat,self.i0Hat,self.r0Hat,self.d0Hat = (self.s0/self.N), (self.i0/self.N), (self.r0/self.N), 0
        self.tHat = self.ts*self.beta

    def simulate(self,vacc = False):
        if vacc:
            deriv = lambda y,t: self.nondimDerivVacc(y,t,self.rho,self.alpha)
        else:
            deriv = lambda y,t: self.nondimDeriv(y,t,self.rho,self.alpha)
        state0 = self.s0Hat,self.i0Hat,self.r0Hat,self.d0Hat
        sirSim = odeint(deriv,state0,self.ts)
        self.S, self.I, self.R, self.D = sirSim.T
        return self.ts, self.S, self.I, self.R, self.D

    def nondimDeriv(self,y, t, rho,alpha):
        S, I, R,D = y
        dSdt = - S * I
        dIdt = S * I - rho * I - alpha*I
        dRdt = rho * I
        dDdt = alpha*I
        return dSdt, dIdt, dRdt, dDdt

    def nondimDerivVacc(self,y, t, rho,alpha):
        S, I, R,D = y
        dSdt = - S * I - np.heaviside(I-0.001,0)*0.01*S
    #    dSdt = - S * I - 0.01*S
        dIdt = S * I - rho * I - alpha*I
        dRdt = rho * I + np.heaviside(I-0.001,0)*0.01*S
        dDdt = alpha*I
        return dSdt, dIdt, dRdt, dDdt

