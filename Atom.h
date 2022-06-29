#pragma once

#include <array>

class Atom
{
public:
	Atom();
	Atom(double x, double y, double z);
	Atom(double x, double y, double z, double mass);
	~Atom();

	void setPos(double x, double y, double z) { m_position = { x, y, z }; }
	void setPos(std::array<double, 3>& pos) { m_position = pos; }
	std::array<double, 3>& getPos() { return m_position; }

	void setOldForce(double fx, double fy, double fz) { m_old_force = { fx, fy, fz }; }
	void resetOldForce() { m_old_force = { 0.0f, 0.0f, 0.0f }; }
	std::array<double, 3>& getOldForce() { return m_old_force; }
	void makeForceOld() { m_old_force = m_force; }

	void setForce(double fx, double fy, double fz) { m_force = { fx, fy, fz }; }
	void resetForce() { m_force = { 0.0f, 0.0f, 0.0f }; }
	std::array<double, 3>& getForce() { return m_force; }

	void setVelocity(std::array<double, 3>& vel) { m_velocity = vel; }
	std::array<double, 3>& getVelocity() { return m_velocity; }

	void setMass(double dMass) { m_mass = dMass; }
	double getMass() { return m_mass; }

	//void updatePosition();
	//void updateVelocities();
	//void correctPosition();

	double getKineticE();

	//static void setTimeStep(double dTimeStep) { s_timestep = dTimeStep; }
	//static void setBoxSize(double dBoxSize) { s_boxsize = dBoxSize; }

private:
	double m_mass { -1.0f };
	std::array<double, 3> m_position { 0.0f, 0.0f, 0.0f };
	std::array<double, 3> m_velocity { 0.0f, 0.0f, 0.0f };
	std::array<double, 3> m_old_force { 0.0f, 0.0f, 0.0f };
	std::array<double, 3> m_force { 0.0f, 0.0f, 0.0f };

	//// Simulation parameters, made static here
	//static double s_timestep;
	//static double s_boxsize;
};