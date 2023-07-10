"""
Calculations for a discrete boson,
coupled to a lossy boson that is classically driven
"""


import numpy as np

from property import Properties


# Constants
hbar_eV = 6.582119569e-16		# [eV*s]
hbar = 1						# unitless; thus time has units of [rad/eV]
# Parameters
we = 0.9326 / hbar  			# resonant frequency [rad/s]
wc = we - (20e-3) / hbar
game = (67.5e-3) / hbar			# dissipative rate [rad/s]
gamc = 6.075e-10 * game
g = 0.332e-3 / hbar				# coupling constant [rad/s]
# Fano parameters
eta_e = g ** 2 / ( np.absolute((we - wc) - (1j*game/2)) ** 2 )
Omegac = wc + eta_e*(wc - we)
Gamc = gamc + eta_e*game

# Compute g2
# Define time domain and scale
npts = 4*24000
period_we = 2*np.pi / we			# period of emitter frequency [rad/eV]
period_wc = 2*np.pi / wc
# tstep = period_we					# smallest time interval
tau_c = 2*np.pi / gamc				# cavity coherence time
tstep = period_wc
endpt = npts * tstep
t = np.linspace(0, endpt, npts)
tau = t

# g2, _, _ = Properties(Omegac+(0*Gamc), np.array([0]), np.array([tau[0]])).computeg2()
# print(g2)
# g2, _, _ = Properties(Omegac+(0*Gamc), t, np.array([tau[0]])).computeg2()
# g2 = g2[:,0]  # for tau = 0
g2, _, _ = Properties(Omegac+(0*Gamc), t, np.array([0]), na=1).computeg2_equaltime()
# output
with open("g2.txt", "w", encoding="utf-8") as f:
	f.write('g2 dynamics at tau=0, vacuum initial state, w=Ωc\n')
    # f.write('g2 dynamics at tau=0, w=Ωc+4gamc\n')
	f.write('{0:12s} {1:12s}\n'.format('t [ps]', 'g2 [1]'))
	# f.write('{0:12s} {1:12s}\n'.format('t [tau_c]', 'g2 [1]'))
	for t_idx, time in enumerate(t):
		x = time*hbar_eV*1e+12; y = g2[t_idx].real
		f.write('{0:<12.6f} {1:<12.6f}\n'.format(x, y))
		# f.write('{0:<12.6f} {1:<12.6f}\n'.format(time/tau_c, g2[t_idx].real))



## eighty five characters in a line
#123456789012345678901234567890123456789012345678901234567890123456789012345678901234
