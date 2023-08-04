#ifndef StochasticBeamingEmission_HH
#define StochasticBeamingEmission_HH

#include "PhotonBeamingEmission.hh"

class StochasticBeamingEmission: public PhotonBeamingEmission
{
public:
    StochasticBeamingEmission(EMField* field, double dt, double sampleFrac = 1,
        double eMin = 0, bool track = false);

    virtual ~StochasticBeamingEmission();

    void Interact(Particle *part, ParticleList *partList) const override;
    
};
#endif
