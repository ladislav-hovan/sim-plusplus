#include "calc.h"

double calibrateLJ(const ParamsLJ &sParams)
{
	double dRatio = std::pow(sParams.r_m / sParams.cutoff, 6);

	return sParams.epsilon * dRatio * (dRatio - 2);
}

double calculateLJ(double dDist, const ParamsLJ &sParams)
{
	if (dDist >= sParams.cutoff)
		return 0.0;

	static double s_dOffset = calibrateLJ(sParams);

	double dRatio = std::pow(sParams.r_m / dDist, 6);

	return sParams.epsilon * dRatio * (dRatio - 2) - s_dOffset;
}

double getPeriodicDist(const Atom &cFirst, const Atom &cSecond, double dBoxSize)
{
	double dSquareSum = 0.0f;

	const auto &aFirstPos = cFirst.getPos();
	const auto &aSecondPos = cSecond.getPos();

	for (int nCoord = 0; nCoord < 3; ++nCoord)
	{
		double dDiff = std::fabs(aFirstPos[nCoord] - aSecondPos[nCoord]);
		dSquareSum += ((dDiff < dBoxSize / 2) ? (dDiff * dDiff) : (dBoxSize - dDiff) * (dBoxSize - dDiff));
	}

	return std::sqrt(dSquareSum);
}

double getSignedDiff(const Atom &cFirst, const Atom &cSecond, double dBoxSize, int nCoord)
{
	double dDiff = cSecond.getPos()[nCoord] - cFirst.getPos()[nCoord];

	if (std::fabs(dDiff) > dBoxSize / 2)
		dDiff = ((dDiff > 0) ? dDiff - dBoxSize : dDiff + dBoxSize);

	return dDiff;
}

std::array<double, 3> getTotalMomentum(const vector<Atom> &vAtoms)
{
	std::array<double, 3> a_dMomentum = { 0.0f, 0.0f, 0.0f };

	for (const Atom &atom : vAtoms)
		for (int nCoord = 0; nCoord < 3; ++nCoord)
			a_dMomentum[nCoord] += atom.getVelocity()[nCoord] * atom.getMass();

	return a_dMomentum;
}

double sumPairwise(const vector<double> &vdValues, int nFirst, int nLast)
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

double sumKahan(const vector<double> &vdValues)
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

double sumSimple(const vector<double> &vdValues)
{
	double dSum = 0.0f;
	for (double dValue : vdValues)
		dSum += dValue;

	return dSum;
}