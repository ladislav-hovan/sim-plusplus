#pragma once

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include "Atom.h"
#include "calc.h"

using std::vector;
using std::string;
using vector2d = vector< vector<double> >;

class Output
{
public:
	Output() = delete;
	Output(const string &strEnergyFile, const string &strPositionFile);
	~Output();

	void logEnergies(double dKinetic, double dPotential, bool bPrint = false);
	void logPositions(const vector<Atom> &vAtoms);

private:
	std::ofstream m_EnergyOutput;
	std::ofstream m_PositionOutput;
};