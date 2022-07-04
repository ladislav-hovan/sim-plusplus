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

int main()
{
	// Simulation parameters
	// TODO: Implement loading values from a parameter file of some sort
	InputParams sInput;
	
	// Check the parameters for consistency
	if (sInput.dBoxSize <= 2 * sInput.lj_par.cutoff)
	{
		std::cerr << "The box size needs to be bigger than twice the LJ cutoff distance\n";
		exit (1);
	}

	// Creating simulation object
	Simulation Simulation;

	// Place atoms randomly in the initial box
	// TODO: Implement loading positions from a file (PDB maybe)
	Simulation.generateRandomPositions();

	// Generate initial velocities
	// TODO: Add an option to load velocities from a file
	Simulation.generateVelocities();

	// Remove center of mass motion
	Simulation.removeTranslation(true);

	// Set initial simulation time
	double dTime = 0.0f;  // ps

	// Open log files
	std::ofstream Energies(sInput.strEnergyFile);
	std::ofstream Positions(sInput.strPositionFile);

	// Start the timer for program run
	clock_t cTime = std::clock();

	// Main integration loop (Velocity Verlet)
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

	// Close the log files
	Energies.close();
	Positions.close();

	// Report on time spent
	cTime = std::clock() - cTime;
	std::cout << "\nTotal run time of the integration loop: " << double(cTime) / CLOCKS_PER_SEC << "s\n";

	return 0;
}