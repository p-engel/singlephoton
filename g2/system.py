from math import sqrt
from qutip.qobj import Qobj
from qutip.tensor import tensor
from qutip.operators import destroy, qeye
from qutip.states import basis


class MasterEquation:
	"""
	Defines the hamiltonian and collapse operators of 
	two coupled bosons
	"""
	#constants
	def __init__(self, w):
		"""
		known attributes of the parameters of the system in eV
		"""
		# self.hbar = 6.582119569e-16				# [ev*s]
		self.hbar = 1								# unitless; t has units [1/eV]
		self.we = 0.9326 / self.hbar  				#resonant freq. [rad/s]
		self.wc = self.we - (20e-3) / self.hbar
		self.game = (67.5e-3) / self.hbar			# dissipative rate [rad/s]
		self.gamc = 6.075e-10 * self.game
		self.g = 0.332e-3 / self.hbar				# coupling constant [rad/s]
		self.v = 1*0.002355 * 2*self.game					# drive amplitude [eV/hbar]
		self.w = w / self.hbar						# drive frequency [rad/s]
		self.N = 9					 				# levels in Hilbert space

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


class InitialState:
	"""The initial state of cavity and emitter"""
	def __init__(self, na, nb):
		self.na = na		# initial cavity fock-state number
		self.nb = nb		# initial emitter fock-state number
		self.N = 9			# levels in Hilbert space

	def rho(self):
		"""define initial state of cavity and emitter"""
		rho0 = tensor(basis(self.N, self.na), basis(self.N, self.nb))

		return rho0


## eighty five characters in a line
#123456789012345678901234567890123456789012345678901234567890123456789012345678901234
