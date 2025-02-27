#include "RunManager.hh"

#include "GaussianEMField.hh"
#include "FocusingField.hh"
#include "ContinuousEmission.hh"
#include "StochasticEmission.hh"
#include "NonLinearBreitWheeler.hh"
#include "LorentzPusher.hh"
#include "LandauPusher.hh"
#include "ModifiedLandauPusher.hh"
#include <pybind11/stl.h>

#ifdef USEOPENMP
    #include <omp.h>
#endif

RunManager::RunManager():
m_timeStep(0), m_timeEnd(0), m_field(nullptr), m_pusher(nullptr), 
m_generator(nullptr), m_fieldSet(false), m_physSet(false), m_genSet(false),
m_sampleFrac(1), m_useBW(false)
{
    // Update this if I ever add another unit system
    m_units = new UnitsSystem("SI");
}

RunManager::~RunManager()
{
    delete m_units;
    delete m_field;
    delete m_generator;
    delete m_pusher;
    m_processList.clear();
}

void RunManager::setTime(double timeStep, double timeEnd)
{
    m_timeStep = timeStep / m_units->RefTime();
    m_timeEnd = timeEnd / m_units->RefTime();
}

void RunManager::setField(const std::string& fieldType, double maxField, 
    double wavelength, double duration, double waist, double polarisation,
    const ThreeVector& start, const ThreeVector& focus)
{
    if (fieldType == "gaussian" || fieldType == "Gaussian" || 
        fieldType == "focusing" || fieldType == "Focusing") 
    {
        m_fieldSet = true;
        m_fieldType = fieldType;
    }
    else
    {
        m_fieldSet = false;
        std::cerr << "Error: Unknown field type." << std::endl;
        return;
    }
    m_maxField = maxField / m_units->RefEField();
    m_wavelength = wavelength / m_units->RefLength();
    m_fieldDuration = duration / m_units->RefTime();
    m_waist = waist / m_units->RefLength();
    m_polarisation = polarisation;
    m_start = start / m_units->RefLength();
    m_focus = focus / m_units->RefLength();
}

void RunManager::setPhysics(const std::string& physics)
{
    if (physics == "Classical" || physics == "clasical" ||
        physics == "Semiclassical" || physics == "semiclassical" ||
        physics == "Quantum" || physics == "quantum")
    {
        m_physSet = true;
        m_physics = physics;
    } else
    {
        m_physSet = false;
        std::cerr << "Error: Unkown physics type. Choices are: \"Classical\", " 
         "\"Semiclassical\" or \"Quantum\";" << std::endl;
        return;
    }
}

void RunManager::setGenerator(const std::string& particleType,
    const std::string& energyDist, double energyParam1, double energyParam2,
    double radius, double duration, double divergence, 
    const ThreeVector& position, const ThreeVector& direction)
{
    m_genSet = true;
    m_particleType = particleType;
    m_energyDist = energyDist;
    m_energyParam1 = energyParam1 / m_units->RefEnergy(); 
    m_energyParam2 = energyParam2 / m_units->RefEnergy();
    m_radius = radius / m_units->RefLength();
    m_particleDuration = duration / m_units->RefTime();
    m_divergence = divergence;
    m_position = position / m_units->RefLength();
    m_direction = direction;
}

void RunManager::setSampleFraction(double sampleFrac)
{
    if (sampleFrac >= 0 && sampleFrac <= 1)
    {
        m_sampleFrac = sampleFrac;
    } else
    {
        std::cerr << "Error: Sampling fraction must be between 0 and 1."
                  << std::endl;
    }
}


void RunManager::usePairProduction(bool useBW)
{
    m_useBW = useBW;
}

void RunManager::beamOn(int events, int threads)
{

    // If threads set then check for openMP install
#ifndef USEOPENMP
    if (threads != 1)
    {
        std::cerr << "Error: OpenMP version of QEDCascPy not built. Set thread"
                     " to 1 or build OpenMP."
                  << std::endl;
        return;
    }
#endif

    delete m_pusher;
    delete m_generator;
    delete m_field;
    m_processList.clear();

    // Check for time / field properties
    if (m_timeStep == 0 || !m_fieldSet || !m_physSet || !m_genSet)
    {
        std::cerr << "Error in beamOn(): RunManager.setTime(...), "
                      "RunManager.setGenerator(...) and RunManager.setField(...) "
                      "must all be called before running simulation." 
                   << std::endl;
        return;
    }

    // Set the field
    if (m_fieldType == "gaussian" || m_fieldType == "Gaussian")
    {
        m_field = new GaussianEMField(m_maxField, m_wavelength, m_fieldDuration,
            m_waist, m_polarisation, m_start, m_focus);
    } else
    {
        m_field = new FocusingField(m_maxField, m_wavelength,  m_fieldDuration,
            m_waist, m_polarisation, m_start, m_focus);
    }
    
    // Set the physics
    if (m_physics == "Classical" || m_physics == "classical")
    {
        m_pusher = new LandauPusher(m_field, m_timeStep);
        ContinuousEmission* emission = new ContinuousEmission(m_field,
            m_timeStep, true, m_sampleFrac, false, 0);
        m_processList.push_back(emission);
    } else if (m_physics == "Semiclassical" || m_physics == "semiclassical")
    {
        m_pusher = new ModifiedLandauPusher(m_field, m_timeStep);
        ContinuousEmission* emission = new ContinuousEmission(m_field,
            m_timeStep, false, m_sampleFrac, false, 0);
        m_processList.push_back(emission);
    } else
    {
        m_pusher = new LorentzPusher(m_field, m_timeStep);
        StochasticEmission* emission = new StochasticEmission(m_field,
            m_timeStep, m_sampleFrac, false, 0);
        m_processList.push_back(emission);
    }

    if (m_useBW == true)
    {
        NonLinearBreitWheeler* breitWheeler = new NonLinearBreitWheeler(m_field, 
            m_timeStep, false);
        m_processList.push_back(breitWheeler);
    }


    // set the generator
    m_generator = new SourceGenerator(m_particleType, m_energyDist, events, 
        m_energyParam1, m_energyParam2, m_radius, m_particleDuration, 
        m_divergence, m_position, m_direction);

    m_input_P_X = std::vector<std::vector<double>>(threads);
    m_electron_P_X = std::vector<std::vector<double>>(threads);
    m_positron_P_X = std::vector<std::vector<double>>(threads);
    m_photon_P_X = std::vector<std::vector<double>>(threads);

    long unsigned int nEvents = m_generator->GetSourceNumber();
    
#ifdef USEOPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for
#endif
    for (unsigned int i = 0; i < nEvents; i++) //loop primes
    {
#ifdef USEOPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        ParticleList* event = m_generator->GenerateList();
        // Store inital particle properties
        for (unsigned int j = 0; j < event->GetNPart(); j++) // Loop particles
        {
            m_input_P_X[tid].push_back(event->GetParticle(j)->GetMomentum()[0]
                * m_units->RefMomentum());
            m_input_P_X[tid].push_back(event->GetParticle(j)->GetMomentum()[1]
                * m_units->RefMomentum());
            m_input_P_X[tid].push_back(event->GetParticle(j)->GetMomentum()[2]
                * m_units->RefMomentum());
            m_input_P_X[tid].push_back(event->GetParticle(j)->GetPosition()[0]
                * m_units->RefLength());
            m_input_P_X[tid].push_back(event->GetParticle(j)->GetPosition()[1]
                * m_units->RefLength());
            m_input_P_X[tid].push_back(event->GetParticle(j)->GetPosition()[2]
                * m_units->RefLength());
        }

        double time(0);
        while(time < m_timeEnd) //loop time
        {
            // Push particles and interact
            for (unsigned int j = 0; j < event->GetNPart(); j++) // Loop particles
            {
                m_pusher->PushParticle(event->GetParticle(j));
                for (unsigned int proc = 0; proc < m_processList.size(); proc++) // loop processes
                {
                    m_processList[proc]->Interact(event->GetParticle(j), event);
                }
            }
            time += m_timeStep;
        }
        
        // Store final particle properties
        for (long unsigned int j = 0; j < event->GetNPart(); ++j)
        {
            if (event->GetParticle(j)->GetCharge() == -1)
            {
                m_electron_P_X[tid].push_back(event->GetParticle(j)->GetMomentum()[0] 
                    * m_units->RefMomentum());
                m_electron_P_X[tid].push_back(event->GetParticle(j)->GetMomentum()[1]
                    * m_units->RefMomentum());
                m_electron_P_X[tid].push_back(event->GetParticle(j)->GetMomentum()[2]
                    * m_units->RefMomentum());
                m_electron_P_X[tid].push_back(event->GetParticle(j)->GetPosition()[0]
                    * m_units->RefLength());
                m_electron_P_X[tid].push_back(event->GetParticle(j)->GetPosition()[1]
                    * m_units->RefLength());
                m_electron_P_X[tid].push_back(event->GetParticle(j)->GetPosition()[2]
                    * m_units->RefLength());
                m_electron_P_X[tid].push_back(event->GetParticle(j)->GetWeight());
            } else if (event->GetParticle(j)->GetCharge() == 0)
            {
                m_photon_P_X[tid].push_back(event->GetParticle(j)->GetMomentum()[0]
                    * m_units->RefMomentum());
                m_photon_P_X[tid].push_back(event->GetParticle(j)->GetMomentum()[1]
                    * m_units->RefMomentum());
                m_photon_P_X[tid].push_back(event->GetParticle(j)->GetMomentum()[2]
                    * m_units->RefMomentum());
                m_photon_P_X[tid].push_back(event->GetParticle(j)->GetPosition()[0]
                    * m_units->RefLength());
                m_photon_P_X[tid].push_back(event->GetParticle(j)->GetPosition()[1]
                    * m_units->RefLength());
                m_photon_P_X[tid].push_back(event->GetParticle(j)->GetPosition()[2]
                    * m_units->RefLength());
                m_photon_P_X[tid].push_back(event->GetParticle(j)->GetWeight());
            } else if (event->GetParticle(j)->GetCharge() == 1)
            {
                m_positron_P_X[tid].push_back(event->GetParticle(j)->GetMomentum()[0]
                    * m_units->RefMomentum());
                m_positron_P_X[tid].push_back(event->GetParticle(j)->GetMomentum()[1]
                    * m_units->RefMomentum());
                m_positron_P_X[tid].push_back(event->GetParticle(j)->GetMomentum()[2]
                    * m_units->RefMomentum());
                m_positron_P_X[tid].push_back(event->GetParticle(j)->GetPosition()[0]
                    * m_units->RefLength());
                m_positron_P_X[tid].push_back(event->GetParticle(j)->GetPosition()[1]
                    * m_units->RefLength());
                m_positron_P_X[tid].push_back(event->GetParticle(j)->GetPosition()[2]
                    * m_units->RefLength());
                m_positron_P_X[tid].push_back(event->GetParticle(j)->GetWeight());
            }
        }
        m_generator->FreeSources(event);
    }
}

py::array_t<double> RunManager::getInput()
{
    std::vector<double> accum;
    for (auto& sub : m_input_P_X)
    {
        accum.insert(accum.end(), sub.begin(), sub.end());
    }
    int particles = accum.size() / 6;
    py::array_t<double> result = py::cast(accum);
    result.resize({particles, 6});
    return result;
}

py::array_t<double> RunManager::getElectrons()
{
    std::vector<double> accum;
    for (auto& sub : m_electron_P_X)
    {
        accum.insert(accum.end(), sub.begin(), sub.end());
    }
    int particles = accum.size() / 7;
    py::array_t<double> result = py::cast(accum);
    result.resize({particles, 7});
    return result;
}

py::array_t<double> RunManager::getPositrons()
{
    std::vector<double> accum;
    for (auto& sub : m_positron_P_X)
    {
        accum.insert(accum.end(), sub.begin(), sub.end());
    }
    int particles = accum.size() / 7;
    py::array_t<double> result = py::cast(accum);
    result.resize({particles, 7});
    return result;
}

py::array_t<double> RunManager::getPhotons()
{
    std::vector<double> accum;
    for (auto& sub : m_photon_P_X)
    {
        accum.insert(accum.end(), sub.begin(), sub.end());
    }
    int particles = accum.size() / 7;
    py::array_t<double> result = py::cast(accum);
    result.resize({particles, 7});
    return result;
}
