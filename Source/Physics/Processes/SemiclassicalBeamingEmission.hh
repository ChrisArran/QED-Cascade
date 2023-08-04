#ifndef SemiclassicalBeamingEmission_HH
#define SemiclassicalBeamingEmission_HH

#include "PhotonBeamingEmission.hh"

class SemiclassicalBeamingEmission: public PhotonBeamingEmission
{
public:
    SemiclassicalBeamingEmission(EMField* field, double dt,
        double sampleFrac = 1, double eMin = 0, bool track = false);
    
    virtual ~SemiclassicalBeamingEmission();

    void Interact(Particle* part, ParticleList *partList) const override;
};
#endif
