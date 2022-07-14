#pragma once

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include "Atom.h"
#include "calc.h"

using std::vector;
using std::string;

class Output
{
public:
	Output() = delete;
	Output(const string &strEnergyFile, const string &strPositionFile);
	~Output();

	void logEnergies(vector<Atom> &vAtoms, double dBoxSize, ParamsLJ& sParams, bool bPrint = false);
	void logPositions(vector<Atom> &vAtoms);

private:
	std::ofstream m_EnergyOutput;
	std::ofstream m_PositionOutput;
};