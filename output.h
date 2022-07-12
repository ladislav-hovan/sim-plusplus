#pragma once

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include "Atom.h"
#include "calc.h"

using std::vector;
using std::string;

void logEnergies(vector<Atom> &vAtoms, std::ofstream &Output, double dBoxSize,
	ParamsLJ &sParams, bool bPrint = false);
void logPositions(vector<Atom> &vAtoms, std::ofstream &Output);

class Output
{
public:
	Output() = delete;
	Output(string strEnergyFile, string strPositionFile);
	~Output();

	void logEnergies(vector<Atom>& vAtoms, double dBoxSize, ParamsLJ& sParams, bool bPrint = false);
	void logPositions(vector<Atom>& vAtoms);

private:
	std::ofstream m_EnergyOutput;
	std::ofstream m_PositionOutput;
};