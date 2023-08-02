"""A set of classes for calculating the properties of the system"""
import numpy as np

# from qutip.rhs_generate import rhs_clear
# from qutip.cy.utilities import _cython_build_cleanup
# from qutip.states import ket2dm
from qutip.qobj import Qobj #, isket
from qutip.correlation import correlation_3op_1t
from qutip.correlation import correlation_3op_2t
from qutip.correlation import _correlation_me_2t
from qutip.mesolve import mesolve
from qutip.solver import Options, Result

from system import MasterEquation, InitialState


# class SecondOrderCorrelation:
class Properties:
	""" Define parameters, system, and state"""
	def __init__(self, w, t, tau, na=0, nb=0):
		self.w = w								# drive frequency [eV]
		self.t = t								# time [s]
		self.tau = tau							# delay times >= 0 [s]
		self.na = na							# initial cavity fock-state number
		self.nb = nb							# initial emitter fock-state number
		self.lindblad = MasterEquation(self.w)
		self.state = InitialState(self.na, self.nb)

	def denomg2(self):
		"""
		Compute the normalization factor for
		the non-steadystate second-order correlation function of
		the light transmitted through the discrete mode
		i.e normalize with n(t)*n(t+tau) to obtain g2
		"""
		H, c_ops, a, _ = self.lindblad.operators()
		rho0 = self.state.rho()
		norm_mat = np.zeros([np.size(self.t), np.size(self.tau)], dtype=complex)
		# first calculate the occupation number as a function of initial time
		n_t1 = mesolve(H, rho0, self.t, c_ops, [a.dag()*a]).expect[0]

		for t_idx, tdelay in enumerate(self.tau):
			# next calculate the occupation number as a function of delay time
			n_t2 = mesolve(H, rho0, (self.t + tdelay), c_ops, [a.dag()*a]).expect[0]
			norm_mat[:,t_idx] = n_t1 * n_t2

		# np.outer(n_t, n_tau, out=norm_mat)

		return norm_mat
		
	def computeg2(self):
		"""
		Compute the normalized second-order correlation 
		of light transmitted through the discrete mode
		"""
		H, c_ops, a, _ = self.lindblad.operators()
		rho0 = self.state.rho()
		# rho0 = None
		G2 = correlation_3op_2t(H, rho0, self.t, self.tau,
			c_ops, a.dag(), a.dag()*a, a, solver="me", 
			# options=Options(ntraj=[20, 100]))
			)
		norm_const = self.denomg2()
		g2 = G2 / norm_const

		return g2, G2, norm_const

	def computeg2_equaltime(self, options=Options(ntraj=[20, 100])):
	    """
	    Hacking qutip's internal function "_correlation_me_2t()" for 
	    calculating the three-operator, two-time,
	    correlation function at tau=0:
	    <A(t)B(t)C(t)>
	    using a master equation solver.
	    The solver only works for positive time differences and 
	    the correlator require positive tau
	    """
	    H, c_ops, a, b = self.lindblad.operators()
	    rho0 = self.state.rho()

	    na = mesolve(H, rho0, self.t, c_ops, [a.dag()*a]).expect[0]
	    nb = mesolve(H, rho0, self.t, c_ops, [b.dag()*b]).expect[0]

	    G2_tau0_a = _correlation_me_2t(H, rho0, self.t, [0], 
			c_ops, a.dag(), a.dag()*a, a)
	    g2_tau0_a = G2_tau0_a[:,0] / (na*na)
	    G2_tau0_b = _correlation_me_2t(H, rho0, self.t, [0], 
			c_ops, b.dag(), b.dag()*b, b)
	    g2_tau0_b = G2_tau0_b[:,0] / (nb*nb)

	    return g2_tau0_a, g2_tau0_b, na, nb


## eighty five characters in a line
#123456789012345678901234567890123456789012345678901234567890123456789012345678901234
