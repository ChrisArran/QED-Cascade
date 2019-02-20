#include "Photon.hh"

Photon::Photon(const ThreeVector &position, const ThreeVector &momentum,
			   double time, bool tracking):
Particle(0.0, 0.0, position, momentum, time, tracking)
{
}

Photon::Photon(double energy, const ThreeVector &position, const ThreeVector &direction,
			   double time, bool tracking):
Particle(0.0, 0.0, position, energy * direction.Norm(), time, tracking)
{
}

Photon::~Photon()
{

}

double Photon::GetEnergy() const
{
	return m_momentum.Mag();
}

double Photon::GetGamma() const
{
	return m_momentum.Mag();
}

double Photon::GetBeta() const
{
	return 1.0;
}

