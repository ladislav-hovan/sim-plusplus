#pragma once

#include <array>

using std::array;

class Atom
{
public:
	Atom();
	Atom(double x, double y, double z);
	Atom(double x, double y, double z, double mass);
	~Atom();

	void setPos(double x, double y, double z) { m_position = { x, y, z }; }
	void setPos(std::array<double, 3>& pos) { m_position = pos; }
	array<double, 3> &getPos() { return m_position; }
	const array<double, 3> &getPos() const { return m_position; }

	void setOldForce(double fx, double fy, double fz) { m_old_force = { fx, fy, fz }; }
	void setOldForce(std::array<double, 3>  &old_force) { m_old_force = old_force; }
	void resetOldForce() { m_old_force = { 0.0f, 0.0f, 0.0f }; }
	array<double, 3>& getOldForce() { return m_old_force; }
	void makeForceOld() { m_old_force = m_force; }

	void setForce(double fx, double fy, double fz) { m_force = { fx, fy, fz }; }
	void resetForce() { m_force = { 0.0f, 0.0f, 0.0f }; }
	array<double, 3>& getForce() { return m_force; }

	void setVelocity(std::array<double, 3> &vel) { m_velocity = vel; }
	array<double, 3>& getVelocity() { return m_velocity; }
	const array<double, 3>& getVelocity() const { return m_velocity; }

	void setMass(double dMass) { m_mass = dMass; }
	double getMass() const { return m_mass; }

	double getKineticE() const;

private:
	double m_mass { -1.0f };
	array<double, 3> m_position { 0.0f, 0.0f, 0.0f };
	array<double, 3> m_velocity { 0.0f, 0.0f, 0.0f };
	array<double, 3> m_old_force { 0.0f, 0.0f, 0.0f };
	array<double, 3> m_force { 0.0f, 0.0f, 0.0f };
};