#pragma once

#include <fstream>
#include <iostream>
#include <iomanip>

#include "Atom.h"
#include "calc.h"

using std::vector;

void logEnergies(vector<Atom> &vAtoms, std::ofstream &Output, double dBoxSize,
	ParamsLJ &sParams, bool bPrint = false);
void logPositions(vector<Atom> &vAtoms, std::ofstream &Output);