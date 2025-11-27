#include <cmath>
#include <fstream>

#include "PhotonBeamingEmission.hh"
#include "Photon.hh"
#include "Numerics.hh"
#include "MCTools.hh"
#include "UnitsSystem.hh"


PhotonBeamingEmission::PhotonBeamingEmission(EMField* field, double dt, double sampleFrac, 
    double eMin, bool track): PhotonEmission(field, dt,sampleFrac, eMin, track)
{
    LoadTables();
    LoadAngleTables();
}

double PhotonBeamingEmission::CalculateZ(double eta, double u) const
{
    // More complex, 2D interpolation in eta and u, using a parallelogram:
    // (etaIndex,uIndexlow)->(etaIndex,uIndexlow+1) : lowlow->lowhigh
    // (etaIndex+1,uIndexhigh)->(etaIndex+1,uIndexhigh+1) : highlow->highhigh
    double rand = MCTools::RandDouble(0, 1);
    int etaIndex, uIndexlow, uIndexhigh;
    double fracEta, fracUlow, fracUhigh;
    double lowValue, highValue, logz, z;
    if (std::log10(eta) > m_phZ_etaAxis[0])
    {
        try {
            Numerics::ClosestPoints(m_phZ_etaAxis, m_phZ_etaLength, std::log10(eta),
                etaIndex, fracEta);
        catch (int queryPoint) {
            std::cerr << "Error: Failed in PhotonBeamingEmission::CalculateZ, with log10(eta)="
                << std::log10(eta) << std::endl;
            std::exit(-1);
        }
        }
        if (std::log10(u) > m_phZ_uAxis[etaIndex][0])
        {
            try {
                Numerics::ClosestPoints(m_phZ_uAxis[etaIndex], m_phZ_uLength, std::log10(u),
                    uIndexlow, fracUlow);
                Numerics::ClosestPoints(m_phZ_uAxis[etaIndex+1], m_phZ_uLength, std::log10(u),
                    uIndexhigh, fracUhigh);
            catch (int queryPoint) {
                std::cerr << "Error: Failed in PhotonBeamingEmission::CalculateZ, with log10(u)="
                    << std::log10(u) << std::endl;
                std::exit(-1);                
            }
            double lowlowValue = Numerics::Interpolate1D(m_phZ_dataTable[etaIndex][uIndexlow],
                m_phZ_zAxis[etaIndex][uIndexlow], m_phZ_zLength, rand);
            double lowhighValue = Numerics::Interpolate1D(m_phZ_dataTable[etaIndex][uIndexlow+1],
                m_phZ_zAxis[etaIndex][uIndexlow+1], m_phZ_zLength, rand);
            double highlowValue = Numerics::Interpolate1D(m_phZ_dataTable[etaIndex+1][uIndexhigh],
                m_phZ_zAxis[etaIndex+1][uIndexhigh], m_phZ_zLength, rand);
            double highhighValue = Numerics::Interpolate1D(m_phZ_dataTable[etaIndex+1][uIndexhigh+1],
                m_phZ_zAxis[etaIndex+1][uIndexhigh+1], m_phZ_zLength, rand);
            lowValue = (1.0 - fracUlow)*lowlowValue + fracUlow*lowhighValue;
            highValue = (1.0 - fracUhigh)*highlowValue + fracUhigh*highhighValue;
            logz = (1.0 - fracEta)*lowValue + fracEta*highValue;
            z = std::pow(10.0, logz);
        }
        else // u < umin
        {
            // Extrapolate downwards in u using a constant u*z/eta
            lowValue = Numerics::Interpolate1D(m_phZ_dataTable[etaIndex][0],
                m_phZ_zAxis[etaIndex][0], m_phZ_zLength, rand);
            highValue = Numerics::Interpolate1D(m_phZ_dataTable[etaIndex+1][0],
                m_phZ_zAxis[etaIndex+1][0], m_phZ_zLength, rand);
            logz = (1.0 - fracEta)*lowValue + fracEta*highValue;
            z = std::pow(10.0,logz + m_phZ_uAxis[etaIndex][0]) / u;
        }
    } else // eta < etamin
    {
        if (std::log10(u) > m_phZ_uAxis[0][0])
        {
            // Extrapolate downwards in eta using a constant u*z/eta
            try {
                Numerics::ClosestPoints(m_phZ_uAxis[0], m_phZ_uLength, std::log10(u),
                        uIndexlow, fracUlow);
            } catch (int queryPoint) {
                std::cerr << "Error: Failed in PhotonBeamingEmission::CalculateZ, with log10(u)="
                    << std::log10(u) << std::endl;
                std::exit(-1);                   
            }
            lowValue = Numerics::Interpolate1D(m_phZ_dataTable[0][uIndexlow],
                    m_phZ_zAxis[0][uIndexlow], m_phZ_zLength, rand);
            highValue = Numerics::Interpolate1D(m_phZ_dataTable[0][uIndexlow+1],
                    m_phZ_zAxis[0][uIndexlow+1], m_phZ_zLength, rand);
            logz = (1.0 - fracUlow)*lowValue + fracUlow*highValue;
            z = std::pow(10.0,logz - m_phZ_etaAxis[0]) * eta;
        }
        else // both eta < etamin and u < umin
        {
            // Extrapolate downwards in both eta and u using a constant u*z/eta
            logz = Numerics::Interpolate1D(m_phZ_dataTable[0][0],
                    m_phZ_zAxis[0][0], m_phZ_zLength, rand);
            z = std::pow(10.0,logz + m_phZ_uAxis[0][0] - m_phZ_etaAxis[0]) * eta / u;
        }
    }
    return z;
}

PhotonBeamingEmission::~PhotonBeamingEmission()
{
    UnloadTables();
    UnloadAngleTables();
}

void PhotonBeamingEmission::LoadAngleTables()
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
    std::ifstream phZFile(path + "/z_sokolov.table");
    std::ifstream phUminmaxFile(path + "/uminmax.table");
    std::ifstream phZminmaxFile(path + "/zminmax.table");

    if (!phZFile)
    {
        std::cerr << "ERROR: Data table for photon angle sampling not found!" 
                  << std::endl;
        std::cerr << "No file at: " << path + "/z_sokolov.table" << std::endl;
        exit(1); 
    }
    if (!phUminmaxFile)
    {
        std::cerr << "ERROR: Data table for photon angle energy grid not found!" 
                  << std::endl;
        std::cerr << "No file at: " << path + "/uminmax.table" << std::endl;
        exit(1); 
    }
    if (!phZminmaxFile)
    {
        std::cerr << "ERROR: Data table for photon angle grid not found!" 
                  << std::endl;
        std::cerr << "No file at: " << path + "/zminmax.table" << std::endl;
        exit(1); 
    }

    // set up photon angle tables
    double logMaxEta, logMinEta;
    phZFile >> m_phZ_etaLength >> m_phZ_uLength >> m_phZ_zLength >> logMinEta >> logMaxEta;
    m_phZ_etaAxis = new double [m_phZ_etaLength];
    m_phZ_uAxis = new double* [m_phZ_etaLength];
    m_phZ_zAxis = new double** [m_phZ_etaLength];
    m_phZ_dataTable = new double** [m_phZ_etaLength];
    for (unsigned int i = 0; i < m_phZ_etaLength; i++)
    {
        m_phZ_uAxis[i] = new double [m_phZ_uLength];
        m_phZ_zAxis[i] = new double* [m_phZ_uLength];
        m_phZ_dataTable[i] = new double* [m_phZ_uLength];
        for (unsigned int j = 0; j < m_phZ_uLength; j++)
        {
            m_phZ_zAxis[i][j] = new double [m_phZ_zLength];
            m_phZ_dataTable[i][j] = new double [m_phZ_zLength];
        }
    }

    // load data
    double deltaEta = (logMaxEta - logMinEta) / (m_phZ_etaLength - 1.0);
    double logMinU,logMaxU,deltaU;
    double logMinZ,logMaxZ,deltaZ,dummyZ;
    for (unsigned int i = 0; i < m_phZ_etaLength; i++)
    {
        m_phZ_etaAxis[i] = logMinEta + i * deltaEta;
        phUminmaxFile >> logMinU >> logMaxU;
        deltaU = (logMaxU - logMinU) / (m_phZ_uLength - 1.0);

        for (unsigned int j = 0; j < m_phZ_uLength; j++)
        {
            m_phZ_uAxis[i][j] = logMinU + j * deltaU;
            phZminmaxFile >> logMinZ >> logMaxZ;
            deltaZ = (logMaxZ - logMinZ) / (m_phZ_zLength - 1.0);

            for (unsigned int k = 0; k < m_phZ_zLength; k++)
            {
                m_phZ_zAxis[i][j][k] = logMinZ + k * deltaZ;
                phZFile >> m_phZ_dataTable[i][j][k];
            }
        }
    }

    phZFile.close();
    phUminmaxFile.close();
    phZminmaxFile.close();

}

void PhotonBeamingEmission::UnloadAngleTables()
{
    for (unsigned int i = 0; i < m_phEn_etaLength; i++)
    {
        for (unsigned int j = 0; j < m_phEn_etaLength; j++)
        {
            delete [] m_phZ_zAxis[i][j];
            delete [] m_phZ_dataTable[i][j];
        }
        delete [] m_phZ_uAxis[i];
        delete [] m_phZ_zAxis[i];
        delete [] m_phZ_dataTable[i];
    }
    delete [] m_phZ_etaAxis;
    delete [] m_phZ_uAxis;
    delete [] m_phZ_zAxis;
    delete [] m_phZ_dataTable;
}
