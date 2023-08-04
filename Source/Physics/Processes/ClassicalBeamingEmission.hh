#ifndef ClassicalBeamingEmission_HH
#define ClassicalBeamingEmission_HH

#include "PhotonEmission.hh"

class ClassicalBeamingEmission: public PhotonEmission
{
public:
    ClassicalBeamingEmission(EMField* field, double dt,
        double sampleFrac = 1, double eMin = 0, bool track = false);
    
    virtual ~ClassicalBeamingEmission();

    void Interact(Particle* part, ParticleList *partList) const override;

private:
    void LoadTables(); // Overload these from PhotonEmission and simplify load-in

    void UnloadTables();

    double CalculateZ() const;

    double CalculateU(double eta, double z) const;

    double* m_phZeta_dataTable;
    double* m_phZeta_zetaAxis;
    unsigned int m_phZeta_zetaLength;
};
#endif
