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
gamc = 6.075e-20 * game
g = 0.332e-3 / hbar				# coupling constant [rad/s]
# Fano parameters
eta_e = g ** 2 / ( np.absolute((we - wc) - (1j*game/2)) ** 2 )
Omegac = wc + eta_e*(wc - we)
Gamc = gamc + eta_e*game
q = (Omegac - wc + (1j*gamc/2)) / (Gamc/2)
epsq_plus = (1 - np.absolute(q)**2) / (2*q.real) \
    + np.sqrt(((1 - np.absolute(q)**2) / (2*q.real))**2 + 1)
epsq_minus = (1 - np.absolute(q)**2) / (2*q.real) \
    - np.sqrt(((1 - np.absolute(q)**2) / (2*q.real))**2 + 1)

# important resonant frequencies across Fano profile
# for negative q
wpeak = epsq_minus*(Gamc/2) + Omegac
wdip = epsq_plus*(Gamc/2) + Omegac
wres = Omegac
wlist = np.array([wpeak, wres, wdip])

# Compute g2
# Define time domain and scale
# npts = 4*24000
npts = 300
period_we = 2*np.pi / we			# period of emitter frequency [rad/eV]
period_wc = 2*np.pi / wc
tauc = 2*np.pi / gamc				# cavity coherence time
tstep = (1/1) * period_we
endpt = npts * tstep
t = np.linspace(0, endpt, npts)
tau = t

# g2, _, _ = Properties(Omegac+(0*Gamc), t, np.array([tau[0]])).computeg2()
# memory allocation ?
g2names = ['g2_wpeak', 'g2_wres', 'g2_wdip']

for w_idx, w in enumerate(wlist):
    g2a, g2b, na, nb = Properties(w, t, np.array([0]), 
        na=1, nb=1).computeg2_equaltime()
    # output
    filename = '../data/' + g2names[w_idx] + '.txt'
    with open(filename, 'w', encoding='utf-8') as f:
    	f.write('g2 dynamics at tau=0, state:Fock-1\n')
    	f.write('{0:13s} {1:14s} {2:14s} {3:13s} {4:13s}\n'.format('t [ps]',
        'g2a [unitless]', 'g2b [unitless]','na [unitless]', 'nb [unitless]'))
    	for t_idx, time in enumerate(t):
    		time_ps = time*hbar_eV*1e+12
    		f.write('{0:<13.6f} {1:<14.6f} {2:<14.6f} '\
                '{3:<13.6f} {4:<13.6f}\n'.format(time_ps,
            g2a[t_idx].real, g2b[t_idx].real, na[t_idx], nb[t_idx]))


## eighty five characters in a line...it's actually seventy nine
#123456789012345678901234567890123456789012345678901234567890123456789012345678901234
