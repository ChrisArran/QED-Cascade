from QEDCascPy import RunManager, ThreeVector
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.constants import pi,e,m_e,c

### Time variables
d_t = 0.05e-15

### Set Field Variables
field_type = "gaussian"
wavelength = 0.8e-6
polarisation = 0

duration_fwhm = 35.0e-15
waist_fwhm = 2*1.22*0.8e-6

### Set particle source
particle_type = "Electron"
energy_dist = "gaussian"
spot_size = 0.1e-6
duration_part = 0.
direction  = ThreeVector(0, 0, 1)
l0 = 0.0

divergence = 1e-3 / np.sqrt(8.0*np.log(2))

mec2 = m_e*c**2
mev = 1e6*e
mev_c = mev/c

def a0_to_field(a0,wavelength):
        omega = 2*pi*c / wavelength
        E = a0 * m_e*omega*c/e
        return E

mev = 1e6*e


def sim_rr_angle(a0,gamma,gamma_spread,
		l0=0, angle=0,x_offset=0,y_offset=0,z_offset=0,
		duration_fwhm = 35.0e-15, waist_fwhm = 2*1.22*0.8e-6,
		Nsamples=10000,up_scale=1,photon_fraction=1.0, reject = 0.0, model="Quantum"):

	Efield = a0_to_field(a0,wavelength)
	energy_mean = gamma*mec2
	energy_sig = gamma_spread*mec2

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
	    energy_sig, spot_size, duration_part, divergence, position,  direction, l0, reject);

	run_manager.setPhysics(model)
	run_manager.usePairProduction(False, up_scale)
	run_manager.setSampleFraction(photon_fraction)

	### Run simulation
	run_manager.beamOn(events = int(Nsamples), threads = 4)

	### Get inputs / outputs
	initial = run_manager.getInput()
	final = run_manager.getElectrons()
	photons = run_manager.getPhotons()

	initial_E = np.sqrt(initial[:,0]**2+initial[:,1]**2+initial[:,2]**2+(0.511*mev_c)**2) - (0.511*mev_c)
	final_E = np.sqrt(final[:,0]**2+final[:,1]**2+final[:,2]**2+(0.511*mev_c)**2) - (0.511*mev_c)
	photons_E = np.sqrt(photons[:,0]**2+photons[:,1]**2+photons[:,2]**2)
	photons_w = photons[:,6]

	rr_rate = 1 - np.sum(final_E)/np.sum(initial_E)
	photon_rate = np.sum(photons_E*photons_w)/np.sum(initial_E)

	return photon_rate


def sim_positrons_angle(a0,gamma,gamma_spread,z_offset=0,x_offset=0,y_offset=0,l0=0,angle=0,
		duration_fwhm = 35.0e-15, waist_fwhm = 2*1.22*0.8e-6,
		Nsamples=10000,up_scale=1,photon_fraction=1.0, reject=0.0, model="Quantum",
		e_folder='electron_data/', g_folder='photon_data/', p_folder='positron_data/'):

	Efield = a0_to_field(a0,wavelength)
	energy_mean = gamma*mec2
	energy_sig = gamma_spread*mec2

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
	    energy_sig, spot_size, duration_part, divergence, position,  direction, l0, reject)

	run_manager.setPhysics(model)
	run_manager.usePairProduction(True, up_scale)
	run_manager.setSampleFraction(photon_fraction)

	### Run simulation
	run_manager.beamOn(events = int(Nsamples), threads = 4)

	### Get inputs / outputs
#	photons = run_manager.getPhotons()
	positrons = run_manager.getPositrons()

#	photons_E = np.sqrt(photons[:,0]**2+photons[:,1]**2+photons[:,2]**2)/mev_c
#	positron_E = np.sqrt((0.511*mev_c)**2+positrons[:,0]**2+positrons[:,1]**2+positrons[:,2]**2)/mev_c
#	photons_weight = photons[:,6]
	positrons_weight = positrons[:,6]

#	bins_g = np.histogram_bin_edges(photons_E[~np.isnan(photons_E)], bins='doane')
#	hist_g,bins_g = np.histogram(photons_E, weights=photons_weight, bins=bins_g )
#	hist_g = hist_g / (Nsamples*photon_fraction)
#	bins_g = 0.5*(bins_g[0:-1]+bins_g[1:])
	
#	g_file = g_folder+'%0.1f_%0.0f_%0.2f_%0.1f_ghist.txt' % (a0, gamma, gamma_spread/gamma, z_offset*1e6)
#	np.savetxt(g_file,np.column_stack([bins_g,hist_g]),header='# E (MeV)\t Count')

	if len(positrons)>0:
		positron_rate = np.sum(positrons_weight)/Nsamples
#		bins_p = np.histogram_bin_edges(positron_E[~np.isnan(positron_E)], bins='doane')
#		hist_p,bins_p = np.histogram(positron_E, weights=positrons_weight, bins=bins_p)
#		hist_p = hist_p / Nsamples
#		bins_p = 0.5*(bins_p[0:-1]+bins_p[1:])
#
#		p_file = p_folder+'%0.1f_%0.0f_%0.2f_%0.1f_%0.1e_phist.txt' % (a0, gamma, gamma_spread/gamma, z_offset*1e6, up_scale)
#		np.savetxt(p_file,np.column_stack([bins_p,hist_p]), header='# E (MeV)\t Count')
	else:
		positron_rate = 0

	return positron_rate


