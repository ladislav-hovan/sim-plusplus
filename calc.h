#pragma once

#include <random>
#include <cmath>
#include <vector>
#include <iostream>

#include "params.h"
#include "Atom.h"

using std::vector;

const double c_dBoltzmann = 1.38064852e-23;  // J/K
const double c_dMu = 1.660539040e-27;  // kg
const double c_dKToNatural = c_dBoltzmann / (1e6 * c_dMu);

double getRand();

double getPeriodicDist(Atom &cFirst, Atom &cSecond, double dBoxSize);
double calculateLJ(double dDist, ParamsLJ &sParams);
void updateForces(Atom &cFirst, Atom &cSecond, double dBoxSize, ParamsLJ &sParams);

double sumPairwise(vector<double> &vdValues, int nFirst = 0, int nLast = -1);
double sumKahan(vector<double> &vdValues);
double sumSimple(vector<double> &vdValues);

vector<Atom> generateRandomPositions(int nAtoms, double dBoxSize, double dMass, double dLimit = 0.2);
void generateVelocities(vector<Atom> &vAtoms, double dTemp);
std::array<double, 3> getTotalMomentum(vector<Atom> &vAtoms);
void removeTranslation(vector<Atom> &vAtoms, bool bReport = false);