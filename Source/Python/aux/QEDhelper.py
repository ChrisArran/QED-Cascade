#!\user\bin\python3
# Chris Arran
# September 2021
# Helper function to define constants etc

from scipy.constants import pi,e,m_e,c

def a0_to_field(a0,wavelength):
	omega = 2*pi*c / wavelength
	E = a0 * m_e*omega*c/e
	return E

mev = 1e6*e
