#pragma once

#include <vector>
#include <array>

class Atom
{
public:
	Atom();
	Atom(double x, double y, double z);
	Atom(double x, double y, double z, double mass);
	~Atom();

	void setPos(double x, double y, double z) { m_position = { x, y, z }; }
	std::array<double, 3>& getPos() { return m_position; }

	void setOldForces(double fx, double fy, double fz) { m_old_forces = { fx, fy, fz }; }
	void resetOldForces() { m_old_forces = { 0.0f, 0.0f, 0.0f }; }
	std::array<double, 3>& getOldForces() { return m_old_forces; }

	void setForces(double fx, double fy, double fz) { m_forces = { fx, fy, fz }; }
	void resetForces() { m_forces = { 0.0f, 0.0f, 0.0f }; }
	std::array<double, 3>& getForces() { return m_forces; }

	std::array<double, 3>& getVelocities() { return m_velocities; }

	void setMass(double dMass) { m_mass = dMass; }
	double getMass() { return m_mass; }

	void updatePosition();
	void updateVelocities();
	void correctPosition();

	double getKineticE();

	static void setTimeStep(double dTimeStep) { s_timestep = dTimeStep; }
	static void setBoxSize(double dBoxSize) { s_boxsize = dBoxSize; }

private:
	double m_mass = -1.0;
	std::array<double, 3> m_position { 0.0f, 0.0f, 0.0f };
	std::array<double, 3> m_velocities { 0.0f, 0.0f, 0.0f };
	std::array<double, 3> m_old_forces{ 0.0f, 0.0f, 0.0f };
	std::array<double, 3> m_forces { 0.0f, 0.0f, 0.0f };

	// Simulation parameters, made static here
	static double s_timestep;
	static double s_boxsize;
};