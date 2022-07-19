#include "main.h"

int main(int argc, char* argv[])
{
	// TODO: Create an Input object, with the ability to check input parameters
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

	// Start the timer for program run
	clock_t cTime = std::clock();

	// TODO: Report properly on step zero

	// Main integration loop (Velocity Verlet)
	// The initial time and step are set by the Simulation constructor
	for (; !Simulation.isFinished(); Simulation.advanceTime())
	{
		if (Simulation.getStep() % 1000 == 0)
			std::cout << "Step " << Simulation.getStep() << "\n";

		if (Simulation.getStep() % 100 == 0)
			Simulation.logPositions();

		// Update the position according to Verlet algorithm
		Simulation.updatePositions();
		Simulation.correctPositions();  // Enforce PBC

		// Update forces for all atoms
		Simulation.updateForces();

		// Update the velocities now
		Simulation.updateVelocities();

		// Calculate the kinetic and potential energies of the system, record them
		if (Simulation.getStep() % 10 == 0)
		{
			Simulation.calculateEnergies();
			Simulation.logEnergies(Simulation.getStep() % 1000 == 0);
		}
	}

	// Report on time spent
	cTime = std::clock() - cTime;
	std::cout << "\nTotal run time of the integration loop: " << double(cTime) / CLOCKS_PER_SEC << "s\n";

	return 0;
}