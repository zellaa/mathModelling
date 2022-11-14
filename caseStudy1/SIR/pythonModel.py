##
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
##
# Define disease parameters
beta, gamma = 0.15, 0.04
rho = gamma/beta
deathRate = 30/100000
alpha = deathRate/beta
# Define initial variables, with "hat" variables those being nondimensionalised
N = 10000000
i0, r0 = 100, 0
s0 = N-i0-r0
s0Hat, i0Hat, r0Hat, d0Hat = (s0/N), (i0/N), (r0/N), 0
y0Hat = s0Hat, i0Hat, r0Hat, d0Hat
t = np.linspace(0, 365, 365)
tHat = t*beta
##
# Define 2 ODE systems to compare a model including vaccination to one without
def undimenDeriv(y, t, rho,alpha):
    S, I, R,D = y
    dSdt = - S * I
    dIdt = S * I - rho * I - alpha*I
    dRdt = rho * I
    dDdt = alpha*I
    return dSdt, dIdt, dRdt, dDdt
def undimenDerivWithVacc(y, t, rho,alpha):
    S, I, R,D = y
    dSdt = - S * I - np.heaviside(I-0.001,0)*0.01*S
#    dSdt = - S * I - 0.01*S
    dIdt = S * I - rho * I - alpha*I
    dRdt = rho * I + np.heaviside(I-0.001,0)*0.01*S
    dDdt = alpha*I
    return dSdt, dIdt, dRdt, dDdt

##
# Integrate the SIR equations over the time grid, t.
sirParams = odeint(undimenDeriv, y0Hat, tHat, args=(rho,alpha))
sirParamsVaccine = odeint(undimenDerivWithVacc, y0Hat, tHat, args=(rho,alpha))
S, I, R, D= sirParams.T
SV, IV, RV, DV = sirParamsVaccine.T
##
# Plot initial nondimensional model
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S, 'b', alpha=0.5, lw=2, label='Susceptible People')
ax.plot(t, I, 'r', alpha=0.5, lw=2, label='Infected People')
ax.plot(t, R, 'g', alpha=0.5, lw=2, label='Recovered People')
#ax.plot(t, D, 'k', alpha=0.5, lw=2, label='Dead People')
# Now plot vaccine model
ax.plot(t, SV, 'c', alpha=0.5, lw=2, label='Susceptible People (Vaccine Model)')
ax.plot(t, IV, 'm', alpha=0.5, lw=2, label='Infected People (Vaccine Model)')
ax.plot(t, RV, 'y', alpha=0.5, lw=2, label='Recovered People (Vaccine Model)')
#ax.plot(t, DV, 'navy', alpha=0.5, lw=2, label='Dead People (Vaccine Model)')
# Define plot parameters
ax.set_xlabel('Time (days)')
ax.set_ylabel('Percentage of Population')
ax.set_ylim(0,1)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(visible=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.title(f"SIR diagram for β = {beta}, γ = {gamma} \n Population of {N} with {i0} initial infected")
plt.show()
##
fig2 = plt.figure(facecolor='w')
ax2 = fig2.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax2.plot(t, D, 'k', alpha=0.5, lw=2, label='Dead People')
# Now plot vaccine model
ax2.plot(t, DV, 'navy', alpha=0.5, lw=2, label='Dead People (Vaccine Model)')
# Define plot parameters
ax2.set_xlabel('Time (days)')
ax2.set_ylabel('Percentage of Population')
#ax.set_ylim(0,1)
ax2.yaxis.set_tick_params(length=0)
ax2.xaxis.set_tick_params(length=0)
ax2.grid(visible=True, which='major', c='w', lw=2, ls='-')
legend = ax2.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.title(f"SIR diagram for β = {beta}, γ = {gamma} \n Population of {N} with {i0} initial infected")
plt.show()
