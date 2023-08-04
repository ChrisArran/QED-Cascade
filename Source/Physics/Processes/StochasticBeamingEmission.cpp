#include "StochasticBeamingEmission.hh"
#include <Numerics.hh>
#include "Photon.hh"
#include <MCTools.hh>

StochasticBeamingEmission::StochasticBeamingEmission(EMField* field, double dt,
    double sampleFrac, double eMin, bool track):
PhotonBeamingEmission(field, dt, sampleFrac, eMin, track)
{
}

StochasticBeamingEmission::~StochasticBeamingEmission()
{
}

void StochasticBeamingEmission::Interact(Particle *part, ParticleList *partList) const
{
    if (part->GetMass() == 0 || part->IsAlive() == false) return;
    // First we need to update the optical depth of the particle based on local values
    double eta = CalculateEta(part);

    // Check for very small values of eta and skip interpolation
    double logh;
    if (std::log10(eta) < m_h_etaAxis[0])
    {
        logh = m_h_dataTable[0];
    } else
    {
        logh = Numerics::Interpolate1D(m_h_etaAxis, m_h_dataTable,
            m_h_length, std::log10(eta));
    }
    double deltaOD = m_dt * std::sqrt(3) * UnitsSystem::alpha * eta
        * std::pow(10.0, logh)
        / (part->GetGamma() * 2.0 * UnitsSystem::pi);
    part->UpdateOpticalDepth(deltaOD);
    // Now check if process hass occured. If so then emmit and react
    if (part->GetOpticalDepth() < 0.0)
    {
        // Calculate photon energy
        double chi = CalculateChi(eta);
        double gammaE = 2.0 * chi * part->GetGamma() / eta;

        // Calculate photon angle
        double u = gammaE / (part->GetGamma() - gammaE);
        double z = CalculateZ(eta, u);
        double betacostheta = 1.0 - std::pow(z,2.0/3.0)/(2.0*std::pow(part->GetGamma(),2.0));
        double theta = acos( betacostheta/std::sqrt(1.0 - std::pow(part->GetGamma(),-2.0)) );
        double phi = MCTools::RandDouble(0, 2.0 * UnitsSystem::pi);
        
        ThreeVector znorm = ThreeVector(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
        ThreeMatrix rot = znorm.RotateToAxis(ThreeVector(0,0,1));
        ThreeVector direction = rot * part->GetDirection();

        ThreeVector gammaP = gammaE * direction;
        part->UpdateTrack(part->GetPosition(), part->GetMomentum() - gammaP);

        // Add new partles to the simulation
        if (gammaE > m_eMin && MCTools::RandDouble(0, 1) < m_sampleFrac)
        {
            Photon* photon = new Photon(gammaE, part->GetPosition(), 
                direction, part->GetWeight() / m_sampleFrac,
                part->GetTime(), m_track);
            partList->AddParticle(photon);
        }
        part->InitOpticalDepth();
    }
}
