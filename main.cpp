#include <iostream>
#include <vector>
#include <ctime>
#include <fstream>

#include "params.h"
#include "calc.h"
#include "Atom.h"
#include "Simulation.h"
#include "output.h"

using std::vector;

//double Atom::s_timestep = 0.0f;
//double Atom::s_boxsize = 0.0f;

int main()
{
	// SImulation parameters
	// TODO: Implement loading values from a parameter file of some sort
	InputParams sInput;
	
	//// LJ parameters, default values used for testing
	//ParamsLJ sParams;

	//// Initial conditions
	//
	//// TODO: Add support for non-square boxes
	//double dBoxSize = 10.0f;  // Length of the side of the simulation box, in nm
	//int nAtoms = 100;
	//double dMass = 39.9623831225f;  // Atomic mass units (Argon-40)
	//double dTemp = 20.0f;  // K
	//sParams.epsilon *= c_dKToNatural;

	if (sInput.dBoxSize <= 2 * sInput.lj_par.cutoff)
	{
		std::cerr << "The box size needs to be bigger than twice the LJ cutoff distance\n";
		exit (1);
	}

	//// Run parameters
	//int nSteps = 100000;
	//double dTimeStep = 0.001f;  // ps

	// Creating simulation object
	Simulation Simulation;

	//// Passing the necessary parameters to atoms
	//Atom::setTimeStep(dTimeStep);
	//Atom::setBoxSize(dBoxSize);

	// Place atoms randomly in the initial box
	// TODO: Implement loading positions from a file (PDB maybe)
	// vector<Atom> vAtoms = generateRandomPositions(nAtoms, dBoxSize, dMass);
	Simulation.generateRandomPositions();

	// Generate initial velocities
	// TODO: Add an option to load velocities from a file
	//generateVelocities(vAtoms, dTemp);
	Simulation.generateVelocities();

	// Remove center of mass motion
	//removeTranslation(vAtoms, true);
	Simulation.removeTranslation(true);

	// Main integration loop (Velocity Verlet)
	double dTime = 0.0f;  // ps
	std::ofstream Energies("energies.dat");
	std::ofstream Positions("traj.dat");
	clock_t cTime = std::clock();
	for (int nStep = 1; nStep <= Simulation.getMaxSteps(); ++nStep, dTime += Simulation.getTimeStep())
	{
		if (nStep % 1000 == 0)
			std::cout << "Step " << nStep << "\n";

		if (nStep % 100 == 0)
			logPositions(Simulation.getAtoms(), Positions);

		// Update the position according to Verlet algorithm
		Simulation.updatePositions();
		Simulation.correctPositions();  // Enforce PBC

		// TODO: Make part of the simulation object
		// Loop over all pairs of atoms and update forces
		for (int nFirst = 0; nFirst < sInput.nAtoms; ++nFirst)
		{
			for (int nSecond = nFirst + 1; nSecond < sInput.nAtoms; ++nSecond)
				updateForces(Simulation.getAtoms()[nFirst], Simulation.getAtoms()[nSecond], sInput.dBoxSize, sInput.lj_par);
		}

		// Update the velocities now
		Simulation.updateVelocities();

		// Calculate the kinetic and potential energies of the system, record them
		if (nStep % 10 == 0)
			logEnergies(Simulation.getAtoms(), Energies, sInput.dBoxSize, sInput.lj_par, (nStep % 1000 == 0));
	}
	Energies.close();
	Positions.close();

	// Report on time spent
	cTime = std::clock() - cTime;
	std::cout << "\nTotal run time of the integration loop: " << double(cTime) / CLOCKS_PER_SEC << "s\n";

	return 0;
}