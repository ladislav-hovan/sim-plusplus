#include "output.h"

void logEnergies(vector<Atom> &vAtoms, std::ofstream &Output, double dBoxSize, 
	ParamsLJ &sParams, bool bPrint)
{
	Output << std::fixed << std::setprecision(4);

	vector<double> vdKineticE;
	vdKineticE.reserve(vAtoms.size());
	for (Atom &atom : vAtoms)
		vdKineticE.push_back(atom.getKineticE());
	double dKinetic = sumPairwise(vdKineticE);

	vector<double> vdPotentialE;
	vdPotentialE.reserve(vAtoms.size());
	for (unsigned int nFirst = 0; nFirst < vAtoms.size(); ++nFirst)
	{
		for (unsigned int nSecond = nFirst + 1; nSecond < vAtoms.size(); ++nSecond)
			vdPotentialE.push_back(calculateLJ(getPeriodicDist(vAtoms[nFirst], vAtoms[nSecond], dBoxSize),
				sParams));
	}
	double dPotential = sumPairwise(vdPotentialE);
	Output << dKinetic << " " << dPotential << " " << dKinetic + dPotential << "\n";

	if (bPrint)
	{
		std::cout << std::fixed << std::setprecision(4);
		std::cout << dKinetic << " " << dPotential << " " << dKinetic + dPotential << "\n";
	}
}

void logPositions(vector<Atom> &vAtoms, std::ofstream &Output)
{
	Output << std::fixed << std::setprecision(4);
	for (Atom &atom : vAtoms)
	{
		for (double fCoord : atom.getPos())
			Output << fCoord << " ";
		Output << "\n";
	}
	Output << "\n";
}