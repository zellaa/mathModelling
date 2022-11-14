from SIR_simulation import *

sir1 = SIR()
t,S,I,R,D = sir1.simulate()
sir2 = SIR()
t,SV,IV,RV,DV = sir2.simulate(vacc=True)


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
plt.title(f"SIR diagram for β = {sir1.beta}, γ = {sir1.gamma} \n Population of {sir1.N} with {sir1.i0} initial infected")
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
plt.title(f"SIR diagram for β = {sir2.beta}, γ = {sir2.gamma} \n Population of {sir2.N} with {sir2.i0} initial infected")
plt.show()
