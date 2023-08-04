#include <fstream>

#include "ContinuousNonLinearBreitWheeler.hh"
#include "Lepton.hh"
#include "Numerics.hh"
#include "UnitsSystem.hh"
#include "MCTools.hh"

ContinuousNonLinearBreitWheeler::ContinuousNonLinearBreitWheeler(EMField* field, double dt, bool track, double up_scale):
NonLinearBreitWheeler(field, dt, track, up_scale)
{
}

ContinuousNonLinearBreitWheeler::~ContinuousNonLinearBreitWheeler()
{
}

void ContinuousNonLinearBreitWheeler::Interact(Particle *part, ParticleList *partList) const
{
    if (part->GetType() != "Photon" || part->IsAlive() == false) return;

    double chi = CalculateChi(part);
    double logt = Numerics::Interpolate1D(m_t_chiAxis, m_t_dataTable,
        m_t_length, std::log10(chi));
    double deltaOD = m_dt * UnitsSystem::alpha * chi * std::pow(10.0, logt)
        / part->GetEnergy();

    // Use m_up_scale < 1 as a sample fraction, so I limit the number of pairs produced
    // Otherwise m_up_scale > 1 to over-sample, creating several pairs per photon per timestep
    if (MCTools::RandDouble(0, 1) > m_up_scale) return;

    double split, pEnergy, eEnergy;
    ThreeVector pMomentum, eMomentum;
    for (int i=0; i<m_up_scale; i++)
    {
        split = CalculateSplit(chi);
        pEnergy = split * part->GetEnergy();
        eEnergy = (1.0 - split) * part->GetEnergy();
        pMomentum =  std::sqrt(pEnergy * pEnergy - 1.0) * part->GetDirection();
        eMomentum =  std::sqrt(eEnergy * eEnergy - 1.0) * part->GetDirection();
        Lepton* positron = new Lepton(1.0, 1.0, part->GetPosition(), pMomentum,
            part->GetWeight()*deltaOD/m_up_scale, part->GetTime(), m_track); 
        Lepton* electron = new Lepton(1.0, -1.0, part->GetPosition(), eMomentum, 
            part->GetWeight()*deltaOD/m_up_scale, part->GetTime(), m_track);
        partList->AddParticle(positron);
        partList->AddParticle(electron);
    }
    
    if (deltaOD < 1.0)
    {
        double reduction = 1.0 - deltaOD;
        Photon* photon = new Photon(part->GetEnergy(), part->GetPosition(), 
                part->GetDirection(), part->GetWeight()*reduction, part->GetTime(), m_track);
        partList->AddParticle(photon);
    }

    part->Kill();
}

