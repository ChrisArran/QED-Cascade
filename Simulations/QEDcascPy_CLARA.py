from QEDCascPy import RunManager, ThreeVector
import numpy as np
from scipy.constants import pi,e,m_e,c

### Time variables
d_t = 0.05e-15

### Set Field Variables
field_type = "gaussian"
wavelength = 0.8e-6 # Laser wavelength, 800nm
polarisation = 0 # 0 for linear horizontal

duration_fwhm = 23.0e-15
waist_fwhm = 2*1.22*0.8e-6

### Set particle source
particle_type = "Electron"
energy_dist = "gaussian"
spot_size = 20e-6
duration_part = 50e-15
direction  = ThreeVector(0, 0, 1) # Propagating along positive z

l0 = 0.0 # Stand-off distance between accelerator and collision
divergence = 1.0e-3 # Electron beam divergence - standard deviation, in radians

# Constants
MeV_c = 1e6*e/c

def a0_to_field(a0,wavelength):
        omega = 2*pi*c / wavelength
        E = a0 * m_e*omega*c/e
        return E

def sim_particledist_angle(a0,gamma,gamma_spread,z_offset=0,x_offset=0,y_offset=0,l0=0,angle=0,
		duration_fwhm = 35.0e-15, waist_fwhm = 2*1.22*0.8e-6, bins_e = [], bins_g = [], bins_p = [],
		Nsamples=10000,up_scale=1,photon_fraction=1.0, model="Quantum"):

	Efield = a0_to_field(a0,wavelength)
	energy_mean = gamma*m_e*c**2
	energy_sig = gamma_spread*m_e*c**2

	w0 = waist_fwhm / np.sqrt(2*np.log(2))
	tau = duration_fwhm / np.sqrt(2*np.log(2))

	min_sep = 2.2*tau*c # start laser at least this far away from electrons
					   # this cuts off at E/E_max < 1%

	focus = ThreeVector(0, 0, 0)
	sinx = np.sin(angle*np.pi/180)
	cosx = np.cos(angle*np.pi/180)
	displacement = x_offset*sinx + z_offset*cosx

	if displacement > -0.5*min_sep: # default laser start is already upstream of focus
		position = ThreeVector(x_offset,y_offset,z_offset - 0.5*min_sep)
		laser_start = ThreeVector((0.5*min_sep+displacement)*sinx,0,(0.5*min_sep+displacement)*cosx)
		end_time = min_sep / c
	else: # make laser start upstream of focus
		position = ThreeVector(x_offset,y_offset,2*z_offset)
		laser_start = ThreeVector(1e-10*sinx,0,1e-10*cosx)
		end_time = (0.5*min_sep+np.abs(displacement)) / c


	### Set inputs
	run_manager = RunManager()
	run_manager.setTime(d_t, end_time)
	run_manager.setField(field_type, Efield, wavelength, 
	    tau, w0, polarisation, laser_start, focus)
	run_manager.setGenerator(particle_type, energy_dist, energy_mean, 
	    energy_sig, spot_size, duration_part, divergence, position,  direction, l0);

	run_manager.setPhysics(model)
	run_manager.usePairProduction(True, up_scale)
	run_manager.setSampleFraction(photon_fraction)

	### Run simulation
	run_manager.beamOn(events = int(Nsamples), threads = 4)

	### Get inputs / outputs
	electrons = run_manager.getElectrons()
	photons = run_manager.getPhotons()
	positrons = run_manager.getPositrons()

	electrons_E = np.sqrt((0.511*MeV_c)**2+electrons[:,0]**2+electrons[:,1]**2+electrons[:,2]**2)/MeV_c
	electrons_weight = electrons[:,6]
	photons_E = np.sqrt(photons[:,0]**2+photons[:,1]**2+photons[:,2]**2)/MeV_c
	photons_weight = photons[:,6]
	positron_E = np.sqrt((0.511*MeV_c)**2+positrons[:,0]**2+positrons[:,1]**2+positrons[:,2]**2)/MeV_c
	positrons_weight = positrons[:,6]

	if len(bins_e)==0:
		bins_e = np.histogram_bin_edges(electrons_E[~np.isnan(electrons_E)], bins='doane')
	hist_e,bins_e = np.histogram(electrons_E, weights=electrons_weight, bins=bins_e )
	hist_e = hist_e / Nsamples
	bins_e = 0.5*(bins_e[0:-1]+bins_e[1:])

	if len(bins_g)==0:
		bins_g = np.histogram_bin_edges(photons_E[~np.isnan(photons_E)], bins='doane')
	hist_g,bins_g = np.histogram(photons_E, weights=photons_weight, bins=bins_g )
	hist_g = hist_g / (Nsamples*photon_fraction)
	bins_g = 0.5*(bins_g[0:-1]+bins_g[1:])

	if len(positrons)>0:
		positron_rate = np.sum(positrons_weight)/Nsamples
		if len(bins_p)==0:
			bins_p = np.histogram_bin_edges(positron_E[~np.isnan(positron_E)], bins='doane')
		hist_p,bins_p = np.histogram(positron_E, weights=positrons_weight, bins=bins_p)
		hist_p = hist_p / Nsamples
		bins_p = 0.5*(bins_p[0:-1]+bins_p[1:])

	else:
		positron_rate = 0
		if len(bins_p)==0:
			bins_p = np.array([0])
		else:
			bins_p = 0.5*(bins_p[0:-1]+bins_p[1:])
		hist_p = np.zeros_like(bins_p)

	return hist_e,bins_e,hist_g,bins_g,hist_p,bins_p
