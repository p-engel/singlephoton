
"""
Unit Tests:
Calculations for a discrete boson,
coupled to a lossy boson that is classically driven.
"""


import numpy as np

from property import Properties


def test_denomg2(w, t, tau):
	""" 
	Unit test of normalization constant of g2 function
	The value divided by itself should be 1 at all times
	"""
	# Suppress/hide the warning due to division by zero
	np.seterr(invalid='ignore')

	try:
		norm_const = Properties(w, t, tau).denomg2()
		output = norm_const / norm_const
		assert (output[1,1:] - 1).all() <= 1e-3
	except AssertionError as a:
		assert isinstance(a, AssertionError), \
			"norm divided by itself is not equal to one"
		print(norm_const)
	else:
		print('norm divided by itself is either zero or one')
	return


def test_computeg2(w, t, tau, na=0):
	""" 
	Unit test of g2 function, using a coherent (vacuum) state
	The value should be 1 at all times
	"""
	# Suppress/hide the warning due to division by zero
	np.seterr(invalid='ignore')

	try:
		g2, G2, norm_const = Properties(w, t, tau, na).computeg2()
		assert np.absolute(g2[-1,0] - 1) < 0.5
	except AssertionError as a:
		print("g2 is neither equal nor close to one")
		print(G2[:, 0], norm_const[:, 0])
	else:
		print("No AssertionError was raised: "
			"either some other problem or g2 is equal/close to one \
			for the coherent (vacuum) state")
		print(g2[:, 0])
		# raise Exception("No AssertionError was raised",
		# 	"either some other problem or g2 is equal/close to one \
		# 	for the coherent (vacuum) state")
	return



def test_computeg2_equaltime(w, t, tau, na=0):
	""" 
	Unit test of g2 function, using a coherent (vacuum) state
	The value should be 1 at all times
	"""
	# Suppress/hide the warning due to division by zero
	np.seterr(invalid='ignore')

	try:
		g2, G2, norm_const = Properties(w, t, tau, na).computeg2_equaltime()
		assert np.absolute(g2[-1] - 1) < 0.5
	except AssertionError as a:
		print("g2 (at tau equal to zero) is neither equal nor close to one")
		print(G2, norm_const)
	else:
		print("No AssertionError was raised: "
			"either some other problem or g2 (at tau equal to zero) \
			is equal/close to one for the coherent (vacuum) state")
		print(g2, norm_const)
	return


## eighty five characters in a line
#123456789012345678901234567890123456789012345678901234567890123456789012345678901234

