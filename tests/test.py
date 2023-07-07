"""
Unit Test:
Calculations for a discrete boson,
coupled to a lossy boson that is classically driven.
"""


import numpy as np

from property import Properties
from unit_tests import test_computeg2, test_denomg2, test_computeg2_equaltime


# Constants
# hbar = 6.582119569e-16		# [ev*s]
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

# Define time domain and scale
# t = np.linspace( 0, (2*np.pi / gamc)*1e-6, 3000 )
npts = 150
period_we = 2*np.pi / we			# period of emitter frequency [rad/eV]
tstep = period_we					# smallest time interval
endpt = npts * tstep
t = np.linspace(0, endpt, npts)
tau = t

# tests
# test_denomg2(Omegac, t, tau)
test_computeg2(Omegac, np.array([0]), np.array([0]))
# test_computeg2_equaltime(Omegac, t, tau)


## eighty five characters in a line
#123456789012345678901234567890123456789012345678901234567890123456789012345678901234
