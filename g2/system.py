from cmath import pi, sqrt
from numpy import linspace

from qutip.qobj import Qobj
from qutip.tensor import tensor
from qutip.operators import destroy, qeye
from qutip.states import basis


class Parameters():
	"""define the parameters for two coupled oscillators"""
	def __init__(self, w):
		"""
		the parameters of the system in eV, based on 
		the measured Fano system published in Nano Lett. by the Masiello group
		"""
		# self.hbar = 6.582119569e-16				# [ev*s]
		self.hbar = 1								# unitless; t has units [1/eV]
		self.we = 0.9326 / self.hbar  				# resonant freq. [rad/s]
		self.wc = self.we - (20e-3) / self.hbar
		self.game = (67.5e-3) / self.hbar			# dissipative rate [rad/s]
		self.gamc = 6.075e-10 * self.game
		self.g = 0.332e-3 / self.hbar				# coupling constant [rad/s]
		self.v = 1*0.002355 * 2*self.game			# drive amplitude [eV/hbar]
		self.w = w / self.hbar			 				# levels in Hilbert space

	def resonances(self):
		"""define important resonant frequencies across the Fano profile"""
		# Fano parameters
		eta_e = g**2 / ((self.w - we)**2 + (game/2)**2)
		Omegac = wc + eta_e*(self.w - we)
		Gamc = gamc + eta_e*game
		q = (Omegac - (wc - 1j*gamc/2)) / (Gamc/2)
		# solutions to the first derivative of fano function equal to zero
		term = (1 - abs(q)**2) / (2*q.real)
		epsq_plus = term + sqrt(term**2 + 1)
		epsq_minus = term - sqrt(term**2 + 1)
		# important resonant frequencies across Fano profile
		# for negative q
		wpeak = epsq_minus*(Gamc/2) + Omegac
		wdip = epsq_plus*(Gamc/2) + Omegac
		wres = Omegac
		wlist = [wpeak, wres, wdip]

		return wlist

	def times(self):
		"""Define time domain and scale"""
		npts = 300
		T_we = 2*pi / we					# period [rad/eV]
		T_wc = 2*pi / wc
		tau_wc = 2*pi / gamc				# coherence time
		tau_we = 2*pi / game
		tstep = period_we
		endpt = npts * tstep
		t = linspace(0, endpt, npts)
		tau = t

		return t, tau


class MasterEquation(Parameters):
	"""
	Defines the hamiltonian and collapse operators of 
	two coupled bosons
	"""
	#constants
	def __init__(self, w):
		"""
		known attributes of the parameters of the system
		"""
		super.__init__(w)
		self.N = 9					# levels in Hilbert space
		self.na = na				# initial cavity fock-state number
		self.nb = nb				# initial emitter fock-state number

	def operators(self):
		""" 
		define the hamiltonian and collapse operators of 
		two coupled bosons
		"""
		a = tensor(destroy(self.N), qeye(self.N))
		b = tensor(qeye(self.N), destroy(self.N))
		H = self.hbar * (
			(self.wc - self.w)*a.dag()*a
			+ (self.we - self.w)*b.dag()*b
			+ self.g*(a.dag()*b + b.dag()*a)
			+ self.v*(b + b.dag())
			)
		Cops = [sqrt(2*self.gamc)*a, sqrt(2*self.game)*b]

		return H, Cops, a, b

	def rho(self):
		"""define initial state of cavity and emitter"""
		rho0 = tensor(basis(self.N, self.na), basis(self.N, self.nb))

		return rho0



## seventy-nine characters in a line
#123456789012345678901234567890123456789012345678901234567890123456789012345678
