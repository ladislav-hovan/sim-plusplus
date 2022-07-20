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
using array2d = array< array<double, 3>, 3 >;

class Simulation
{
public:
	Simulation() = delete;
	Simulation(const InputParams& sInput);
	~Simulation() = default;

	double getTimeStep() const { return m_dTimeStep; }
	double getBoxSize() const { return m_dBoxSize; }
	int getMaxSteps() const { return m_nMaxSteps; }

	void generateRandomPositions(double dLimit = 0.2);
	void generateVelocities();
	void loadPositions(const string &strPosFile);
	void loadVelocities(const string &strVelFile);
	void removeTranslation(bool bReport = false);

	void updatePositions();
	void updateVelocities();
	void updateForces();
	void correctPositions();

	void logEnergies(bool bPrint = false) { m_Output.logEnergies(m_dKineticE, m_dPotentialE, bPrint); }
	void logPositions() { m_Output.logPositions(m_vAtoms); }

	void advanceTime() { m_nStep += 1; m_dTime += m_dTimeStep; }
	bool isFinished() const { return m_nStep >= m_nMaxSteps; }
	int getStep() const { return m_nStep; }

	void updateDistances();

	void calculateEnergies();
	double getPotentialE() { return m_dPotentialE; }
	double getKineticE() { return m_dKineticE; }
	double getTotalE() { return m_dPotentialE + m_dKineticE; }

	void reportStep() { std::cout << "Step " << m_nStep << "\n"; }

private:
	// Vector of all the atoms
	vector<Atom> m_vAtoms{};
	
	// Vector of distances between all the atoms
	vector2d m_vvdDistances{};

	// Box vectors
	array2d m_aadBoxVectors{};

	// Vectors and doubles to hold potential and kinetic energy
	vector<double> m_vdPotentialE{};
	vector<double> m_vdKineticE{};
	double m_dPotentialE{ 0.0f };
	double m_dKineticE{ 0.0f };

	// A class responsible for file output
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

	// Current state of the simulation
	double m_dTime{ 0.0f };  // ps
	int m_nStep{ 0 };

	// PRNG
	int m_nSeed{ -1 };
	std::mt19937 m_Mersenne{};

	// Private functions, not to be called from outside
	// PRNG related
	void initialisePRNG();
	double getRand();

	// Initialisation of the Simulation object
	void initialiseDistances();
	void copyForcesToOld();
	void resetForces();
	void initialiseDistancesAndForces();
	void initialiseEnergyVectors();

	// Energy calculation
	void calculatePotentialE();
	void calculateKineticE();
};