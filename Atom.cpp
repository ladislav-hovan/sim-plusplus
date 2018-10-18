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

void Atom::updatePosition()
{
	for (int nCoord = 0; nCoord < 3; ++nCoord)
	{
		m_position[nCoord] += m_velocities[nCoord] * s_timestep;
		m_position[nCoord] += m_old_forces[nCoord] * std::pow(s_timestep, 2) / (2 * m_mass);
	}

	// TODO: Add some check to warn against big changes in position
}

void Atom::updateVelocities()
{
	for (int nCoord = 0; nCoord < 3; ++nCoord)
	{
		m_velocities[nCoord] += (m_old_forces[nCoord] + m_forces[nCoord]) * s_timestep / (2 * m_mass);
		m_old_forces[nCoord] = m_forces[nCoord];
	}

	resetForces();
}

void Atom::correctPosition()
{
	for (int nCoord = 0; nCoord < 3; ++nCoord)
	{
		if (m_position[nCoord] < 0.0)
			m_position[nCoord] += s_boxsize;
		if (m_position[nCoord] >= s_boxsize)
			m_position[nCoord] -= s_boxsize;
	}
}

double Atom::getKineticE()
{
	// In natural units
	double dSquaredVels = 0.0;
	for (auto &vel : m_velocities)
		dSquaredVels += std::pow(vel, 2);

	return 0.5f * m_mass * dSquaredVels;
}