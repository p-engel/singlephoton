"""A set of classes for calculating the properties of the system"""
import numpy as np

from qutip.qobj import Qobj #, isket
from qutip.correlation import correlation_3op_1t
from qutip.correlation import correlation_3op_2t
from qutip.correlation import _correlation_me_2t
from qutip.mesolve import mesolve
from qutip.solver import Options, Result

from system import Parameters MasterEquation


# class SecondOrderCorrelation:
class Properties:
	""" Define parameters, system, and state"""
	def __init__(self, w, na=0, nb=0):
		self.w = w								# drive frequency [eV]
		self.na = na							# initial cavity fock-state number
		self.nb = nb							# initial emitter fock-state number
		self.pars = Parameters(self.w)
		self.lindblad = MasterEquation(self.w)

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
	    times, _ = self.pars.times()
	    taus = [0]
	    H, c_ops, a, b, rho0 = self.lindblad.operators()
	    ops = [a, b]
	    g2list = []
	    nlist = []

	    for it_op, op in enumerate(ops):
	    	# average occupation number
	    	n = mesolve(H, rho0, times, c_ops, [op.dag()*op]).expect[0]
	    	# second order correlation
	    	G2 = _correlation_me_2t(H, rho0, times, taus, 
				c_ops, op.dag(), op.dag()*op, op)
	    	# normalized second order correlation
	    	g2 = G2[:,0] / (n**2)
	    	# add to list
	    	nlist.append()
	    	g2list.append()

	    return g2list, nlist


## eighty five characters in a line
#123456789012345678901234567890123456789012345678901234567890123456789012345678901234
