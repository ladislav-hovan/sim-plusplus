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

// TODO: The energy vectors should be permanent, possibly outside of output class (for virial and stuff)
void Output::logEnergies(vector<Atom> &vAtoms, double dBoxSize, ParamsLJ &sParams, bool bPrint)
{
	vector<double> vdKineticE;
	vdKineticE.reserve(vAtoms.size());
	for (Atom& atom : vAtoms)
		vdKineticE.push_back(atom.getKineticE());
	double dKinetic = sumPairwise(vdKineticE);

	vector<double> vdPotentialE;
	vdPotentialE.reserve(vAtoms.size());
	for (unsigned int nFirst = 0; nFirst < vAtoms.size(); ++nFirst)
	{
		for (unsigned int nSecond = nFirst + 1; nSecond < vAtoms.size(); ++nSecond)
			vdPotentialE.push_back(calculateLJ(getPeriodicDist(vAtoms[nFirst], vAtoms[nSecond], dBoxSize), sParams));
	}
	double dPotential = sumPairwise(vdPotentialE);

	m_EnergyOutput << dKinetic << " " << dPotential << " " << dKinetic + dPotential << "\n";

	if (bPrint)
	{
		std::cout << std::fixed << std::setprecision(4);
		std::cout << dKinetic << " " << dPotential << " " << dKinetic + dPotential << "\n";
	}
}

void Output::logPositions(vector<Atom> &vAtoms)
{
	for (Atom& atom : vAtoms)
	{
		for (double fCoord : atom.getPos())
			m_PositionOutput << fCoord << " ";
		m_PositionOutput << "\n";
	}

	m_PositionOutput << "\n";
}
