#include <cmath>
#include <fstream>

#include "ClassicalBeamingEmission.hh"
#include "Photon.hh"
#include "Numerics.hh"
#include "MCTools.hh"
#include "ThreeVector.hh"
#include "UnitsSystem.hh"


ClassicalBeamingEmission::ClassicalBeamingEmission(EMField* field, double dt,
    double sampleFrac, double eMin, bool track):
PhotonEmission(field, dt, sampleFrac, eMin, track)
{
    LoadTables();
}

ClassicalBeamingEmission::~ClassicalBeamingEmission()
{
    UnloadTables();
}

void ClassicalBeamingEmission::Interact(Particle *part, ParticleList *partList) const
{
    // Check for alive lepton
    if (part->GetMass() == 0 || part->IsAlive() == false) return;

    // prevent from occuring every time step
    if (MCTools::RandDouble(0, 1) > m_sampleFrac) return;

    double eta = CalculateEta(part);
    
    // Changes calculation of photon energy and angle to match Tom Blackburn's classical calculation
    // Calculate photon angle
    double z = CalculateZ();
    double betacostheta = 1.0 - std::pow(z,2.0/3.0)/(2.0*std::pow(part->GetGamma(),2.0));
    double theta = acos( betacostheta/std::sqrt(1.0 - std::pow(part->GetGamma(),-2.0)) );
    double phi = MCTools::RandDouble(0, 2.0 * UnitsSystem::pi);
    
    ThreeVector znorm = ThreeVector(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
    ThreeMatrix rot = znorm.RotateToAxis(ThreeVector(0,0,1));
    ThreeVector direction = rot * part->GetDirection();
    
    // Calculate photon energy
    double u = CalculateU(eta, z);
    double gammaE = part->GetGamma()*part->GetMass()*u;

    // Calculate photon weight
    double logh = m_h_dataTable[0];
    double deltaOD = m_dt * std::sqrt(3) * UnitsSystem::alpha * eta
        * std::pow(10.0, logh)
        / (part->GetGamma() * 2.0 * UnitsSystem::pi);
    double weight = part->GetWeight() * deltaOD / m_sampleFrac;

    // Add new partles to the simulation
    if (gammaE > m_eMin)
    {
        Photon* photon = new Photon(gammaE, part->GetPosition(), 
            direction, weight, part->GetTime(), m_track);
        partList->AddParticle(photon);
    }
}

double ClassicalBeamingEmission::CalculateZ() const
{
    double rand = MCTools::RandDouble(0, 1);
    double delta = acos( (-9.0 + 50.0*rand -25.0*std::pow(rand,2.0))/16.0 );
    double z = std::pow( (2.0+4.0*cos(delta/3.0))/(5.0*(1-rand)), 3.0 );
    return z;
}

double ClassicalBeamingEmission::CalculateU(double eta, double z) const
{
    double rand = MCTools::RandDouble(0, 1);
    double logzeta = Numerics::Interpolate1D(m_phZeta_dataTable,
        m_phZeta_zetaAxis, m_phZeta_zetaLength, rand);
    double zeta = std::pow(10.0, logzeta);
    double u = 3.0*eta*zeta/(2.0*z);
    return u;
}

void ClassicalBeamingEmission::LoadTables() // Simplified, this doesn't need h or ksi_sokolov table
{
    char* tablePath(getenv("QED_TABLES_PATH"));
    if (tablePath == NULL)
    {
        std::cout << "Error: Enviromental variable \"QED_TABLES_PATH\" "; 
        std::cout << "is not set!" << std::endl;
        std::cout <<  "Please set QED_TABLES_PATH to point to tables directory."
                  << std::endl;
    }

    std::string path(tablePath);
    std::ifstream phZetaFile(path + "/zeta_sokolov.table");
    if (!phZetaFile)
    {
        std::cerr << "ERROR: Data table for photon energy sampling not found!" 
                  << std::endl;
        std::cerr << "No file at: " << path + "/zeta_sokolov.table" << std::endl;
        exit(1); 
    }

    // classical photon angle table
    phZetaFile >> m_phZeta_zetaLength;
    m_phZeta_dataTable = new double [m_phZeta_zetaLength];
    m_phZeta_zetaAxis   = new double [m_phZeta_zetaLength];
    double  logZeta, zetaCDF;
    for (unsigned int i = 0; i < m_phZeta_zetaLength; i++)
    {
        phZetaFile >> logZeta >> zetaCDF;
        m_phZeta_zetaAxis[i] = logZeta;
        m_phZeta_dataTable[i] = zetaCDF;
    }

}

void ClassicalBeamingEmission::UnloadTables()
{
    delete [] m_phZeta_dataTable;
    delete [] m_phZeta_zetaAxis;
}
