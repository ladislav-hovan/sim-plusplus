#include <iostream>
#include <vector>
#include <ctime>
#include <fstream>

#include "params.h"
#include "calc.h"
#include "Atom.h"
#include "output.h"

using std::vector;

double Atom::s_timestep = 0.0f;
double Atom::s_boxsize = 0.0f;

int main()
{
	// LJ parameters, default values used for testing
	ParamsLJ sParams;

	// Initial conditions
	double dBoxSize = 10.0f;  // Length of the side of the simulation box, in nm
	int nAtoms = 100;
	double dMass = 39.9623831225f;  // Atomic mass units (Argon-40)
	double dTemp = 20.0f;  // K
	sParams.epsilon *= c_dKToNatural;

	if (dBoxSize <= 2 * sParams.cutoff)
	{
		std::cerr << "The box size needs to be bigger than twice the LJ cutoff distance\n";
		std::getchar();
		exit (1);
	}

	// Run parameters
	int nSteps = 100000;
	double dTimeStep = 0.001f;  // ps

	// Passing the necessary parameters to atoms
	Atom::setTimeStep(dTimeStep);
	Atom::setBoxSize(dBoxSize);

	// Place atoms randomly in the initial box
	vector<Atom> vAtoms = generateRandomPositions(nAtoms, dBoxSize, dMass);

	// Generate initial velocities
	generateVelocities(vAtoms, dTemp);

	// Remove center of mass motion
	removeTranslation(vAtoms, true);

	// Main integration loop (Velocity Verlet)
	double dTime = 0.0f;  // ps
	std::ofstream Energies("energies.dat");
	std::ofstream Positions("traj.dat");
	clock_t cTime = std::clock();
	for (int nStep = 1; nStep <= nSteps; ++nStep, dTime += dTimeStep)
	{
		if (nStep % 1000 == 0)
			std::cout << "Step " << nStep << "\n";

		if (nStep % 100 == 0)
			logPositions(vAtoms, Positions);

		// Update the position according to Verlet algorithm
		for (Atom &atom : vAtoms)
		{
			atom.updatePosition();
			atom.correctPosition();
		}

		// Loop over all pairs of atoms and update forces
		for (int nFirst = 0; nFirst < nAtoms; ++nFirst)
		{
			for (int nSecond = nFirst + 1; nSecond < nAtoms; ++nSecond)
				updateForces(vAtoms[nFirst], vAtoms[nSecond], dBoxSize, sParams);
		}

		// Update the velocities now
		for (Atom &atom : vAtoms)
			atom.updateVelocities();

		// Calculate the kinetic and potential energies of the system, record them
		if (nStep % 10 == 0)
			logEnergies(vAtoms, Energies, dBoxSize, sParams, (nStep % 1000 == 0));
	}
	Energies.close();
	Positions.close();

	// Report on time spent
	cTime = std::clock() - cTime;
	std::cout << "\nTotal run time of the integration loop: " << double(cTime) / CLOCKS_PER_SEC << "s\n";

	std::getchar();
	return 0;
}