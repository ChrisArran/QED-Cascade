#ifndef PhotonBeamingEmission_HH
#define PhotonBeamingEmission_HH

#include "PhotonEmission.hh"

class PhotonBeamingEmission: public PhotonEmission
{
public:
    PhotonBeamingEmission(EMField* field, double dt, double sampleFrac, double eMin, bool track);
    
    virtual ~PhotonBeamingEmission();

    virtual void Interact(Particle *part, ParticleList *partList) const = 0;

protected:

    void LoadAngleTables();

    void UnloadAngleTables();

    double CalculateZ(double eta, double u) const;

    double*** m_phZ_dataTable;
    double*** m_phZ_zAxis;
    double** m_phZ_uAxis;
    double* m_phZ_etaAxis;
    unsigned int m_phZ_zLength, m_phZ_uLength, m_phZ_etaLength;
};
#endif
