#include <cmath>
#include <fstream>

#include "SemiclassicalBeamingEmission.hh"
#include "Photon.hh"
#include "Numerics.hh"
#include "MCTools.hh"
#include "UnitsSystem.hh"


SemiclassicalBeamingEmission::SemiclassicalBeamingEmission(EMField* field, double dt,
    double sampleFrac, double eMin, bool track): 
PhotonBeamingEmission(field, dt, sampleFrac, eMin, track)
{
}

SemiclassicalBeamingEmission::~SemiclassicalBeamingEmission()
{
}

void SemiclassicalBeamingEmission::Interact(Particle *part, ParticleList *partList) const
{
    // Check for alive lepton
    if (part->GetMass() == 0 || part->IsAlive() == false) return;

    // prevent from occuring every time step
    if (MCTools::RandDouble(0, 1) > m_sampleFrac) return;

    double eta = CalculateEta(part);

    // Check for very small values of eta / classical
    double logh, chi;
    if (std::log10(eta) < m_h_etaAxis[0])
    {
        logh = m_h_dataTable[0];
        double chiMin = std::exp(m_h_etaAxis[0]);
        // Extrapolate table for classical (stolen from Chris Arran).
        chi = CalculateChi(chiMin) * (eta / chiMin)  * (eta / chiMin);
    } else
    {
        logh = Numerics::Interpolate1D(m_h_etaAxis, m_h_dataTable,
            m_h_length, std::log10(eta));
        chi = CalculateChi(eta);
    }

    // Calculate photon weight
    double deltaOD = m_dt * std::sqrt(3) * UnitsSystem::alpha * eta
        * std::pow(10.0, logh)
        / (part->GetGamma() * 2.0 * UnitsSystem::pi);
    double weight = part->GetWeight() * deltaOD / m_sampleFrac;

    // Calculate photon energy
    double gammaE = 2.0 * chi * part->GetGamma() / eta;
    ThreeVector gammaP = gammaE * part->GetDirection();

    // Calculate photon angle
    double u = gammaE / (part->GetGamma() - gammaE);
    double z = CalculateZ(eta, u);
    double betacostheta = 1.0 - std::pow(z,2.0/3.0)/(2.0*std::pow(part->GetGamma(),2.0));
    double theta = acos( betacostheta/std::sqrt(1.0 - std::pow(part->GetGamma(),-2.0)) );
    double phi = MCTools::RandDouble(0, 2.0 * UnitsSystem::pi);

    ThreeVector znorm = ThreeVector(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
    ThreeMatrix rot = znorm.RotateToAxis(ThreeVector(0,0,1));
    ThreeVector direction = rot * part->GetDirection();

    // Add new partles to the simulation
    if (gammaE > m_eMin)
    {
        Photon* photon = new Photon(gammaE, part->GetPosition(), 
            direction, weight, part->GetTime(), m_track);
        partList->AddParticle(photon);
    }
}
