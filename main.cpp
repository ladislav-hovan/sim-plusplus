#include "main.h"

int main(int argc, char* argv[])
{
	// Simulation parameters
	InputParams sInput;
	// Read from a parameter file if provided on the command line, else use defaults
	if (argc > 1)
		sInput = readParameters(argv[1]);
	
	// Check the parameters for consistency
	if (sInput.dBoxSize <= 2 * sInput.lj_par.cutoff)
	{
		std::cerr << "The box size needs to be bigger than twice the LJ cutoff distance\n";
		exit (1);
	}

	// Creating simulation object
	Simulation Simulation(sInput);

	// Place atoms randomly in the initial box or load initial positions from a file if provided
	if (sInput.strPosInputFile.empty())
		Simulation.generateRandomPositions();
	else
		Simulation.loadPositions(sInput.strPosInputFile);

	// Generate initial velocities or load them from a file if provided
	if (sInput.strVelInputFile.empty())
		Simulation.generateVelocities();
	else
		Simulation.loadVelocities(sInput.strVelInputFile);

	// Remove center of mass motion
	Simulation.removeTranslation(true);

	// Set initial simulation time
	// TODO: Move into simulation
	double dTime = 0.0f;  // ps

	// Open log files
	// TODO: Make the output streams part of a class which will be inside Simulation
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

		// Update forces for all atoms
		Simulation.updateForces();

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