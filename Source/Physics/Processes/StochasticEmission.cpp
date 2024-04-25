#include "StochasticEmission.hh"
#include <Numerics.hh>
#include "Photon.hh"
#include <MCTools.hh>

StochasticEmission::StochasticEmission(EMField* field, double dt,
    double sampleFrac, double eMin, bool track):
PhotonEmission(field, dt, sampleFrac, eMin, track)
{
}

StochasticEmission::~StochasticEmission()
{
}

void StochasticEmission::Interact(Particle *part, ParticleList *partList) const
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
    // Now check if process has occurred. If so then emit and react
    if ((m_sampleFrac <= 1.0) && (part->GetOpticalDepth() < 0.0))
    {
        double chi = CalculateChi(eta);
        double gammaE = 2.0 * chi * part->GetGamma() / eta;
        if (gammaE > part->GetGamma())
        {
            std::cerr << "Photon energy too high in StochasticEmission::Interact, "
                << "E = " << gammaE << " m_ec^2" << " from eta = " << eta << ", chi = " << chi
                << std::endl;
        }
        ThreeVector gammaP = gammaE * part->GetDirection();
        part->UpdateTrack(part->GetPosition(), part->GetMomentum() - gammaP);
        // Add new particles to the simulation
        // Downsampling : only actually create a subset of these photons
        if (gammaE > m_eMin && MCTools::RandDouble(0, 1) < m_sampleFrac)
        {
            Photon* photon = new Photon(gammaE, part->GetPosition(), 
                part->GetDirection(), part->GetWeight() / m_sampleFrac,
                part->GetTime(), m_track);
            partList->AddParticle(photon);
        }
        part->InitOpticalDepth();
    } else if (m_sampleFrac > 1.0)
    // Upscaling : sample at higher than the normal rate of emission and reduce particle weights accordingly
    // Use a different algorithm to ensure emissions are still independent
    {
        double lambda = deltaOD * m_sampleFrac;
        int n = MCTools::RandPoisson(lambda);
        // Emit a photon
        if (n > 0)
        {
            double chi = CalculateChi(eta);
            double gammaE = 2.0 * chi * part->GetGamma() / eta;
            double newweight = part->GetWeight()*n/m_sampleFrac;
            if (gammaE > m_eMin)
            {
                Photon* photon = new Photon(gammaE, part->GetPosition(), 
                        part->GetDirection(), newweight,
                        part->GetTime(), m_track);
                partList->AddParticle(photon);
            }

            part->SetWeight( part->GetWeight() - newweight );
            ThreeVector gammaP = gammaE * part->GetDirection();

            // Create a new daughter lepton that experiences radiation reaction
            // Enforce a minimum weight to stop exponential increase in particle number
            if (newweight > pow(m_sampleFrac,-3))
            {
                Lepton* littlepart = new Lepton(part->GetMass(), part->GetCharge(), part->GetPosition(), 
                        part->GetMomentum() - gammaP, newweight,
                        part->GetTime(), m_track);
                partList->AddParticle(littlepart);
            }
        }
    }
}

