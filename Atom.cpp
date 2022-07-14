#include "Atom.h"
#include "calc.h"

Atom::Atom()
{
}

Atom::~Atom()
{
}

Atom::Atom(double x, double y, double z)
{
	m_position = { x, y, z };
}

Atom::Atom(double x, double y, double z, double mass)
{
	m_position = { x, y, z };
	m_mass = mass;
}

double Atom::getKineticE() const
{
	// In natural units
	double dSquaredVels = 0.0;

	for (auto &vel : m_velocity)
		dSquaredVels += std::pow(vel, 2);

	return 0.5f * m_mass * dSquaredVels;
}