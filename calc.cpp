#include "calc.h"

double getRand()
{
	//static std::random_device Rand;
	//static std::mt19937 Mersenne(Rand());
	static std::mt19937 Mersenne(1992);  // For testing
	static std::uniform_real_distribution<double> Distribution(0.0, 1.0);

	return Distribution(Mersenne);
}

double calibrateLJ(ParamsLJ &sParams)
{
	double dRatio = std::pow(sParams.r_m / sParams.cutoff, 6);

	return sParams.epsilon * dRatio * (dRatio - 2);
}

double calculateLJ(double dDist, ParamsLJ &sParams)
{
	if (dDist >= sParams.cutoff)
		return 0.0;

	static double dOffset = calibrateLJ(sParams);

	double dRatio = std::pow(sParams.r_m / dDist, 6);

	return sParams.epsilon * dRatio * (dRatio - 2) - dOffset;
}

double getPeriodicDist(Atom &cFirst, Atom &cSecond, double dBoxSize)
{
	double dSquareSum = 0.0f;
	for (int nCoord = 0; nCoord < 3; ++nCoord)
	{
		double dDiff = std::fabs(cSecond.getPos()[nCoord] - cFirst.getPos()[nCoord]);
		dSquareSum += ((dDiff < dBoxSize / 2) ? std::pow(dDiff, 2) : std::pow(dBoxSize - dDiff, 2));
	}

	return std::sqrt(dSquareSum);
}

double getSignedDiff(Atom &cFirst, Atom &cSecond, double dBoxSize, int nCoord)
{
	double dDiff = cSecond.getPos()[nCoord] - cFirst.getPos()[nCoord];
	if (std::fabs(dDiff) > dBoxSize / 2)
		dDiff = ((dDiff > 0) ? dDiff - dBoxSize : dDiff + dBoxSize);

	return dDiff;
}

void updateForces(Atom &cFirst, Atom &cSecond, double dBoxSize, ParamsLJ &sParams)
{
	double dDist = getPeriodicDist(cFirst, cSecond, dBoxSize);
	if (dDist > sParams.cutoff)
		return;

	double dRatio = std::pow(sParams.r_m / dDist, 6);
	double dMagnitude = 12 * sParams.epsilon * dRatio * (1 - dRatio) / std::pow(dDist, 2);

	for (int nCoord = 0; nCoord < 3; ++nCoord)
	{
		double dProduct = dMagnitude * getSignedDiff(cFirst, cSecond, dBoxSize, nCoord);
		cFirst.getForces()[nCoord] += dProduct;
		cSecond.getForces()[nCoord] -= dProduct;
	}
}

std::array<double, 3> getTotalMomentum(vector<Atom> &vAtoms)
{
	std::array<double, 3> a_dMomentum = { 0.0f, 0.0f, 0.0f };
	for (Atom &atom : vAtoms)
		for (int nCoord = 0; nCoord < 3; ++nCoord)
			a_dMomentum[nCoord] += atom.getVelocities()[nCoord] * atom.getMass();
	return a_dMomentum;
}

double sumPairwise(vector<double> &vdValues, int nFirst, int nLast)
{
	if (nLast == -1)
		nLast = static_cast<int>(vdValues.size());
	if ((nLast - nFirst) <= 3)
	{
		double dResult = vdValues[nFirst];
		for (int nPos = nFirst + 1; nPos < nLast; ++nPos)
			dResult += vdValues[nPos];
		return dResult;
	}
	else
	{
		int nMiddle = (nLast + nFirst) / 2;
		return sumPairwise(vdValues, nFirst, nMiddle)
			+ sumPairwise(vdValues, nMiddle, nLast);
	}
}

double sumKahan(vector<double> &vdValues)
{
	double dSum = 0.0f;
	double dDiff = 0.0f;
	for (double dValue : vdValues)
	{
		double dToAdd = dValue - dDiff;
		double dTotal = dSum + dToAdd;
		dDiff = (dTotal - dSum) - dToAdd;
		dSum = dTotal;
	}
	return dSum;
}

double sumSimple(vector<double> &vdValues)
{
	double dSum = 0.0f;
	for (double dValue : vdValues)
		dSum += dValue;

	return dSum;
}

void removeTranslation(vector<Atom> &vAtoms, bool bReport)
{
	std::array<double, 3> a_dMomentum = { 0.0f, 0.0f, 0.0f };
	double dTotalMass = 0.0f;
	for (Atom &atom : vAtoms)
	{
		dTotalMass += atom.getMass();
		for (int nCoord = 0; nCoord < 3; ++nCoord)
			a_dMomentum[nCoord] += atom.getVelocities()[nCoord] * atom.getMass();
	}

	if (bReport)
	{
		std::cout << "Initial momentum vector:\n";
		std::cout << a_dMomentum[0] << " " << a_dMomentum[1] << " " << a_dMomentum[2] << "\n";
	}

	for (double &fPart : a_dMomentum)
		fPart /= dTotalMass;
	for (Atom &atom : vAtoms)
		for (int nCoord = 0; nCoord < 3; ++nCoord)
			atom.getVelocities()[nCoord] -= a_dMomentum[nCoord];
	
	if (bReport)
	{
		a_dMomentum = getTotalMomentum(vAtoms);
		std::cout << "\nMomentum vector after CoM motion removal:\n";
		std::cout << a_dMomentum[0] << " " << a_dMomentum[1] << " " << a_dMomentum[2] << "\n\n";
	}
}

vector<Atom> generateRandomPositions(int nAtoms, double dBoxSize, double dMass, double dLimit)
{
	vector<Atom> vAtoms;
	vAtoms.reserve(nAtoms);

	for (int nCount = 0; nCount < nAtoms; ++nCount)
	{
		bool bAccept = false;
		Atom cTemp;
		while (!bAccept)
		{
			cTemp = Atom(dBoxSize * getRand(), dBoxSize * getRand(), dBoxSize * getRand(),
				dMass);
			bAccept = true;
			for (int nCheck = 0; nCheck < nCount; ++nCheck)
				if (getPeriodicDist(vAtoms[nCheck], cTemp, dBoxSize) < dLimit)
					bAccept = false;
		}
		vAtoms.push_back(cTemp);
	}

	return vAtoms;
}

void generateVelocities(vector<Atom> &vAtoms, double dTemp)
{
	// Velocities have a normal distribution with std = sqrt(kT/m) for each component
	//std::random_device Rand;
	//std::mt19937 Mersenne(Rand());
	std::mt19937 Mersenne(1992);  // For testing

	for (Atom &atom : vAtoms)
	{
		std::normal_distribution<double> Normal(0.0, std::sqrt(dTemp * c_dKToNatural / atom.getMass()));
		for (int nCoord = 0; nCoord < 3; ++nCoord)
			atom.getVelocities()[nCoord] = Normal(Mersenne);
	}
}