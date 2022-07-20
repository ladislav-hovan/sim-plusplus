#include "main.h"

int main(int argc, char* argv[])
{
	// Class for dealing with the input
	Input sInput;
	// Read from a parameter file if provided on the command line, else use defaults
	if (argc > 1)
		sInput.readParameters(argv[1]);

	// Check the input parameters for consistency
	sInput.checkValues();
	
	// Creating simulation object
	Simulation Simulation(sInput.getParams());

	// Place atoms randomly in the initial box or load initial positions from a file if provided
	if (sInput.getParams().strPosInputFile.empty())
		Simulation.generateRandomPositions();
	else
		Simulation.loadPositions(sInput.getParams().strPosInputFile);

	// Generate initial velocities or load them from a file if provided
	if (sInput.getParams().strVelInputFile.empty())
		Simulation.generateVelocities();
	else
		Simulation.loadVelocities(sInput.getParams().strVelInputFile);

	// Remove center of mass motion
	Simulation.removeTranslation(true);

	// Start the timer for program run
	clock_t cTime = std::clock();

	// Report on everything at step zero
	Simulation.reportStep();
	Simulation.logPositions();
	Simulation.calculateEnergies();
	Simulation.logEnergies(true);

	// Main integration loop (Velocity Verlet)
	// The initial time and step are set by the Simulation constructor
	do
	{
		Simulation.advanceTime();

		if (Simulation.getStep() % 1000 == 0)
			Simulation.reportStep();

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
	while (!Simulation.isFinished());

	// Report on time spent
	cTime = std::clock() - cTime;
	std::cout << "\nTotal run time of the integration loop: " << double(cTime) / CLOCKS_PER_SEC << "s\n";

	return 0;
}