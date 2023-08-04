#ifndef CONTINUOUSNONLINEARBREITWHEELER_HH
#define CONTINUOUSNONLINEARBREITWHEELER_HH

#include "NonLinearBreitWheeler.hh"

class ContinuousNonLinearBreitWheeler: public NonLinearBreitWheeler 
{
public:
    ContinuousNonLinearBreitWheeler(EMField* field, double dt, bool track = false, double up_scale = 1.0);

    virtual ~ContinuousNonLinearBreitWheeler();

    void Interact(Particle* part, ParticleList *partList) const override;

private:

};
#endif
