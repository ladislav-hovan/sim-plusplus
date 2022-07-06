#pragma once

#include <random>
#include <cmath>
#include <vector>
#include <iostream>

#include "params.h"
#include "Atom.h"

using std::vector;

double getPeriodicDist(Atom &cFirst, Atom &cSecond, double dBoxSize);
double getSignedDiff(Atom& cFirst, Atom& cSecond, double dBoxSize, int nCoord);
double calculateLJ(double dDist, ParamsLJ &sParams);
void updateForces(Atom &cFirst, Atom &cSecond, double dBoxSize, ParamsLJ &sParams);

double sumPairwise(vector<double> &vdValues, int nFirst = 0, int nLast = -1);
double sumKahan(vector<double> &vdValues);
double sumSimple(vector<double> &vdValues);

std::array<double, 3> getTotalMomentum(vector<Atom> &vAtoms);