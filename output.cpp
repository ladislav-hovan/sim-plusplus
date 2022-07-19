#include "output.h"

Output::Output(const string &strEnergyFile, const string &strPositionFile): m_EnergyOutput{ strEnergyFile }, m_PositionOutput{ strPositionFile }
{
	m_EnergyOutput << std::fixed << std::setprecision(4);
	m_PositionOutput << std::fixed << std::setprecision(4);
}

Output::~Output()
{
	m_EnergyOutput.close();
	m_PositionOutput.close();
}

void Output::logEnergies(double dKinetic, double dPotential, bool bPrint)
{
	// Reporting to file and optionally to standard output
	m_EnergyOutput << dKinetic << " " << dPotential << " " << dKinetic + dPotential << "\n";

	if (bPrint)
	{
		std::cout << std::fixed << std::setprecision(4);
		std::cout << dKinetic << " " << dPotential << " " << dKinetic + dPotential << "\n";
	}
}

void Output::logPositions(const vector<Atom> &vAtoms)
{
	for (const Atom &cAtom : vAtoms)
	{
		for (double dCoord : cAtom.getPos())
			m_PositionOutput << dCoord << " ";

		m_PositionOutput << "\n";
	}

	m_PositionOutput << "\n";
}
