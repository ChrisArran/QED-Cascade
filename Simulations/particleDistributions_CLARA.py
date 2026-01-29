import time
t0 = time.time()
from QEDcascPy_CLARA import sim_particledist_angle
import numpy as np
import matplotlib.pyplot as plt

# Sim parameters
Nsims = 100
up_scale = 1e3
Nsamples = 1e3

# Electron parameters
electron_energy = 250e6 # energy in eV
electron_energyspread = 2.5e6 # energy spread in eV
gamma0 = electron_energy/0.511e6
gamma_spread = electron_energyspread/0.511e6

# Laser parameters
Etot = 8.6 # J
duration0 = 23 # FWHM of intensity (fs)
waist0 = 2 # Beam waist w0 (microns)

angle = 15 # degrees
l0 = 0 # Stand-off distance (meters)
f0 = 0 # Focal shift (meters)

filename = '_%1.1fPW_%1.0fMeV_%1.1fum' % (Etot/duration0,gamma0*0.511e6/1e6,f0*1e6)
params = [l0,f0]; # [length,focus shift]

# Random jitter (standard deviation)
t_std = 25e-15 # in seconds
r_std = 10e-6  # in meters
theta_std = 0  # in radians

# Set up the variables etc
I0 = Etot/((duration0*1e-15)*(waist0*1e-6)**2) * (4*np.log(2)/np.pi)**1.5
a0 = np.sqrt(I0/2.1378e22)

EL = np.sqrt(2*I0/(8.85e-12*3e8))
ES = 1.32e18 # V/m
chi_max = gamma0*EL/ES * (1+np.cos(angle*np.pi/180))
print("Running a simulation with chi_max=",chi_max)

# Create all the random offsets
z_offsets = np.random.default_rng().normal(0, 0.5*t_std*3e8, Nsims)
x_offsets = np.random.default_rng().normal(0, r_std, Nsims)
y_offsets = np.random.default_rng().normal(0, r_std, Nsims)
thetax_offsets = np.random.default_rng().normal(0, theta_std, Nsims)
thetay_offsets = np.random.default_rng().normal(0, theta_std, Nsims)

# Create all the empty arrays to put results in
bins_e = np.linspace(0,4000,401)
bins_g = np.linspace(0,2000,41)
bins_p = np.linspace(0,2000,41)
histe_tot = 0
histg_tot = 0
histp_tot = 0
positronrate = np.zeros(Nsims)
photonrate = np.zeros(Nsims)
rrrate = np.zeros(Nsims)

t1 = time.time()
print("Setup in %0.2f s" % (t1-t0))

# Run the MC function repeatedly to optimize laser duration and waist and write to file
for j in range(Nsims):
	hist_e,bins_e1,hist_g,bins_g1,hist_p,bins_p1 = sim_particledist_angle(a0, gamma0, gamma0*0.1,angle=angle, 
					x_offset=f0*np.sin(np.pi/12)*1e-6+x_offsets[j]+l0*thetax_offsets[j], 
					y_offset=y_offsets[j]+l0*thetay_offsets[j], 
					z_offset=f0*np.cos(np.pi/12)*1e-6+z_offsets[j],
					l0=l0,photon_fraction=1,
					duration_fwhm = duration0*1e-15, waist_fwhm = waist0*1e-6, bins_e=bins_e, bins_g=bins_g, bins_p=bins_p,
					Nsamples=Nsamples, up_scale=up_scale, model="Quantum")

    # Accumulate the total histograms of particle energies over many simulated shots
	histe_tot = histe_tot + hist_e/Nsims
	histg_tot = histg_tot + hist_g/Nsims
	histp_tot = histp_tot + hist_p/Nsims

    # Output the rates of positron and photon production on each simulated shot
	positronrate[j] = np.sum(hist_p)/Nsamples
	photonrate[j] = np.sum(hist_g)/Nsamples
	rrrate[j] = np.sum(bins_g*hist_g)/(Nsamples*gamma0*0.511)

	if (Nsims>=100):
		t2 = time.time()
		print("%i of %i sims complete after %0.2fs" % (j+1,Nsims,t2-t1))

# Mean and standard deviation of the positron rate
mean_rrrate = np.sum(rrrate)/Nsims
std_rrrate = np.sqrt(np.sum(rrrate**2)/Nsims - mean_rrrate**2)

if (Nsims<100):
	t2 = time.time()
	print("%i simulations complete after %0.2fs: " % (Nsims,t2-t1) )

exp = np.floor(np.log10(mean_rrrate))
print('length: %1.3e mm; Radiation rate: (%1.3f+-%1.3f)x10^%1.3f' % (l0*1e3,mean_rrrate/10**exp,std_rrrate/10**exp,exp)

# Write the different particle spectra to text files
dEe = bins_e[1]-bins_e[0]
np.savetxt('electronSpectrum' + filename + '.txt', np.transpose(np.vstack((bins_e[1:],histe_tot/dEe*1e-12/1.6e-19))),header='E (MeV)\tdN/dE (/MeV/pC)')

dEg = bins_g[1]-bins_g[0]
np.savetxt('photonSpectrum' + filename + '.txt', np.transpose(np.vstack((bins_g[1:],histg_tot/dEg*1e-12/1.6e-19))),header='E (MeV)\tdN/dE (/MeV/pC)')

dEp = bins_p[1]-bins_p[0]
np.savetxt('positronSpectrum' + filename + '.txt', np.transpose(np.vstack((bins_p[1:],histp_tot/dEp*1e-12/1.6e-19))),header='E (MeV)\tdN/dE (/MeV/pC)')

# Plot the rate of photon production and their spectra

fig,axs = plt.subplots(2,2)

# Cumulative number of shots producing greater than a certain number of pairs
min_r = np.floor(np.log10(np.min(rrrate)*10e-12/1.6e-19))
if min_r < -6:
	min_r = -6
max_r = np.ceil(np.log10(np.max(rrrate)*10e-12/1.6e-19))
if max_r < 0:
	max_r = 0
bins_r = np.logspace(min_r,max_r,10*int(max_r-min_r)+1)
ratehist,binedges = np.histogram(positronrate*10e-12/1.6e-19,bins=bins_r,weights=np.ones_like(positronrate)/Nsims)
axs[0,0].plot(binedges[:-1],100*np.flipud(np.cumsum(np.flipud(ratehist))))
axs[0,0].set_xscale('log')
axs[0,0].set_xlabel('Radiated Energy > $E_γ/E_-$')
axs[0,0].set_ylabel('Cumulative % of shots')

# Positron spectrum
axs[0,1].stairs(histp_tot*10e-12/1.6e-19, edges=bins_p, fill=True)
axs[0,1].set_xlabel('$E_p$ (MeV)')
axs[0,1].set_ylabel('Pairs per 10 pC')

# Cumulative number of shots producting more than a certain number of photons
bins_gr = np.logspace(-4,2,31)
photonratehist,photonbinedges = np.histogram(photonrate,bins=bins_gr,weights=np.ones_like(photonrate)/Nsims)
axs[1,0].plot(photonbinedges[:-1],100*np.flipud(np.cumsum(np.flipud(photonratehist))))
axs[1,0].set_xscale('log')
axs[1,0].set_xlabel('Photons / electron > $N_γ/N_-$')
axs[1,0].set_ylabel('Cumulative % of shots')

# Photon spectrum
axs[1,1].stairs(histg_tot*10e-12/1.6e-19, edges=bins_g, fill=True)
axs[1,1].set_xlabel('$E_γ$ (MeV)')
axs[1,1].set_ylabel('Photons per 10 pC')
axs[1,1].set_yscale('log')

plt.tight_layout()
plt.savefig('particleDistributions_CLARA_' + filename + '.png')
#plt.show()

