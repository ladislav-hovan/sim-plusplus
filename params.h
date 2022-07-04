#pragma once

#include <string>

#include "constants.h"

using std::string;

struct ParamsLJ
{
	// These data correspond to Argon-40
	double epsilon = 119.8f * Constants::KToNatural;  // Lennard-Jones energy parameter, in natural units after conversion
	double r_m = 1.12246f * 0.3405f;  // Lennard-Jones distance at minimum, in nm
	double cutoff = 2.5f * 0.3405f;  // Cutoff distance of potential function, in nm
};

struct InputParams
{
	// LJ parameters, default values are ok
	ParamsLJ lj_par;

	// Initial conditions for box and atoms
	// TODO: Add support for non-square boxes
	double dBoxSize = 10.0f;  // Length of the side of the simulation box, in nm
	int nAtoms = 100;
	double dMass = 39.9623831225f;  // Atomic mass units (Argon-40)

	// Run parameters
	int nSteps = 100000;
	double dTimeStep = 0.001f;  // ps
	double dTemp = 20.0f;  // K, only used for velocity generation right now

	// Random seed
	int nSeed = -1;  // -1 means a random seed will be used

	// Output files
	string strEnergyFile = "energies.dat";
	string strPositionFile = "traj.dat";
};