#pragma once

#include <vector>
#include <random>

#include "Atom.h"
#include "params.h"
#include "calc.h"
#include "input.h"
#include "Output.h"

using std::vector;
using vectorad = vector< array<double, 3> >;
using vector2d = vector< vector<double> >;

class Simulation
{
public:
	Simulation() = delete;
	Simulation(InputParams& sInput);
	~Simulation();

	double getTimeStep() { return m_dTimeStep; }
	double getBoxSize() { return m_dBoxSize; }
	int getMaxSteps() { return m_nMaxSteps; }

	void generateRandomPositions(double dLimit = 0.2);
	void generateVelocities();
	void loadPositions(const string &strPosFile);
	void loadVelocities(const string &strVelFile);
	void removeTranslation(bool bReport = false);

	void updatePositions();
	void updateVelocities();
	void updateForces();
	void correctPositions();

	void logEnergies(bool bPrint = false) { m_Output.logEnergies(m_vAtoms, m_vvdDistances, m_dBoxSize, m_LJPar, bPrint); }
	void logPositions() { m_Output.logPositions(m_vAtoms); }

	void advanceTime() { m_nStep += 1; m_dTime += m_dTimeStep; }
	bool isFinished() { return m_nStep > m_nMaxSteps; }

	void setTime(double dTime) { m_dTime = dTime; }
	void setStep(int nStep) { m_nStep = nStep; }
	int getStep() { return m_nStep; }

	void updateDistances();

private:
	vector<Atom> m_vAtoms{};

	vector2d m_vvdDistances{};

	Output m_Output;

	// Simulation parameters
	double m_dTimeStep{ 0.001f };  // ps
	double m_dBoxSize{ 10.0f };  // nm
	double m_dTemp{ 20.0f };  // K, only used for initial velocity generation (no thermostat yet)
	int m_nMaxSteps{ 10000 };
	int m_nAtoms{ 100 };

	// Atom properties
	double m_dMass{ 39.9623831225f };  // Atomic mass units (Argon-40)
	ParamsLJ m_LJPar{};  // All the Lennard-Jones parameters

	// Current state of simulation
	double m_dTime{ 0.0f };  // ps
	int m_nStep{ 1 };

	// PRNG
	int m_nSeed{ -1 };
	std::mt19937 m_Mersenne{};

	// Private functions, not to be called from outside
	void initialisePRNG();
	double getRand();
	void initialiseDistances();
};