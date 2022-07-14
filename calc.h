#pragma once

#include <random>
#include <cmath>
#include <vector>
#include <iostream>

#include "params.h"
#include "Atom.h"

using std::vector;

double getPeriodicDist(const Atom &cFirst, const Atom &cSecond, double dBoxSize);
double getSignedDiff(const Atom &cFirst, const Atom &cSecond, double dBoxSize, int nCoord);
double calculateLJ(double dDist, const ParamsLJ &sParams);

double sumPairwise(const vector<double> &vdValues, int nFirst = 0, int nLast = -1);
double sumKahan(const vector<double> &vdValues);
double sumSimple(const vector<double> &vdValues);

std::array<double, 3> getTotalMomentum(const vector<Atom> &vAtoms);