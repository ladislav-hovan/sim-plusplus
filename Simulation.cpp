#include "Simulation.h"

Simulation::Simulation(InputParams& sInput): m_dTimeStep{sInput.dTimeStep}, m_dBoxSize{sInput.dBoxSize}, m_nSeed{sInput.nSeed}, m_nMaxSteps{sInput.nSteps}, m_nAtoms{sInput.nAtoms},
											 m_dMass{sInput.dMass}, m_dTemp{sInput.dTemp}, m_LJPar{sInput.lj_par}, m_Output{sInput.strEnergyFile, sInput.strPositionFile},
											 m_nStep{sInput.nInitialStep}, m_dTime{sInput.dInitialTime}
{
	initialisePRNG();
}

Simulation::~Simulation()
{
}

void Simulation::initialisePRNG()
{
	// If the seed is -1, generate a random one
	if (m_nSeed == -1)
	{
		std::random_device Random;
		m_nSeed = Random();
	}

	// Seed the PRNG
	m_Mersenne = std::mt19937(m_nSeed);
}

double Simulation::getRand()
{
	// Static so only defined once
	static std::uniform_real_distribution<double> Distribution(0.0, 1.0);

	return Distribution(m_Mersenne);
}

void Simulation::generateRandomPositions(double dLimit)
{
	m_vAtoms.clear();
	m_vAtoms.reserve(m_nAtoms);

	for (int nCount = 0; nCount < m_nAtoms; ++nCount)
	{
		bool bAccept = false;
		Atom cTemp;
		while (!bAccept)
		{
			cTemp = Atom(m_dBoxSize * getRand(), m_dBoxSize * getRand(), m_dBoxSize * getRand(), m_dMass);
			bAccept = true;
			for (int nCheck = 0; nCheck < nCount; ++nCheck)
				if (getPeriodicDist(m_vAtoms[nCheck], cTemp, m_dBoxSize) < dLimit)
					bAccept = false;
		}
		m_vAtoms.push_back(cTemp);
	}

	initialiseDistances();
}

void Simulation::generateVelocities()
{
	// Velocities have a normal distribution with std = sqrt(kT/m) for each component
	for (Atom& cAtom : m_vAtoms)
	{
		std::normal_distribution<double> Normal(0.0, std::sqrt(m_dTemp * Constants::KToNatural / cAtom.getMass()));
		for (int nCoord = 0; nCoord < 3; ++nCoord)
			 cAtom.getVelocity()[nCoord] = Normal(m_Mersenne);
	}
}

void Simulation::loadPositions(const string& strPosFile)
{
	vectorad vadValues = readValuesFromFile(strPosFile);

	m_nAtoms = vadValues.size();

	m_vAtoms.clear();
	m_vAtoms.reserve(m_nAtoms);

	for (int nCount = 0; nCount < m_nAtoms; ++nCount)
		m_vAtoms.push_back(Atom(vadValues[nCount][0], vadValues[nCount][1], vadValues[nCount][2], m_dMass));

	initialiseDistances();
}

void Simulation::loadVelocities(const string& strVelFile)
{
	vectorad vadValues = readValuesFromFile(strVelFile);

	if (vadValues.size() != m_nAtoms)
	{
		std::cerr << "The length of the list of velocities doesn't match the number of atoms" << std::endl;
		exit(4);
	}

	for (int nCount = 0; nCount < m_nAtoms; ++nCount)
		m_vAtoms[nCount].setVelocity(vadValues[nCount]);
}

void Simulation::removeTranslation(bool bReport)
{
	std::array<double, 3> a_dMomentum = { 0.0f, 0.0f, 0.0f };
	double dTotalMass = 0.0f;
	for (Atom& cAtom : m_vAtoms)
	{
		dTotalMass += cAtom.getMass();
		for (int nCoord = 0; nCoord < 3; ++nCoord)
			a_dMomentum[nCoord] += cAtom.getVelocity()[nCoord] * cAtom.getMass();
	}

	if (bReport)
	{
		std::cout << "Initial momentum vector:\n";
		std::cout << a_dMomentum[0] << " " << a_dMomentum[1] << " " << a_dMomentum[2] << "\n";
	}

	for (double& dPart : a_dMomentum)
		dPart /= dTotalMass;
	for (Atom& cAtom : m_vAtoms)
		for (int nCoord = 0; nCoord < 3; ++nCoord)
			cAtom.getVelocity()[nCoord] -= a_dMomentum[nCoord];

	if (bReport)
	{
		a_dMomentum = getTotalMomentum(m_vAtoms);
		std::cout << "\nMomentum vector after CoM motion removal:\n";
		std::cout << a_dMomentum[0] << " " << a_dMomentum[1] << " " << a_dMomentum[2] << "\n\n";
	}
}

void Simulation::updatePositions()
{
	for (Atom &cAtom : m_vAtoms)
	{
		std::array<double, 3> &adPosition = cAtom.getPos();
		const std::array<double, 3> &adVelocity = cAtom.getVelocity();
		const std::array<double, 3> &adOldForce = cAtom.getOldForce();

		const static double s_dTimeStepSquared = std::pow(m_dTimeStep, 2);
		double dFactor = s_dTimeStepSquared / (2 * cAtom.getMass());

		for (int nCoord = 0; nCoord < 3; ++nCoord)
			adPosition[nCoord] += adVelocity[nCoord] * m_dTimeStep + adOldForce[nCoord] * dFactor;
	}

	updateDistances();
}

void Simulation::updateVelocities()
{
	for (Atom &cAtom : m_vAtoms)
	{
		std::array<double, 3> &adVelocity = cAtom.getVelocity();
		const std::array<double, 3> &adOldForce = cAtom.getOldForce();
		const std::array<double, 3> &adForce = cAtom.getForce();

		double dFactor = m_dTimeStep / (2 * cAtom.getMass());

		for (int nCoord = 0; nCoord < 3; ++nCoord)
			adVelocity[nCoord] += (adOldForce[nCoord] + adForce[nCoord]) * dFactor;

		cAtom.makeForceOld();
		cAtom.resetForce();
	}
}

void Simulation::updateForces()
{
	for (int nFirst = 0; nFirst < m_nAtoms; ++nFirst)
	{
		for (int nSecond = nFirst + 1; nSecond < m_nAtoms; ++nSecond)
		{
			double dDist = m_vvdDistances[nFirst][nSecond];
			if (dDist > m_LJPar.cutoff)
				continue;

			double dRatio = std::pow(m_LJPar.r_m / dDist, 6);
			double dMagnitude = 12 * m_LJPar.epsilon * dRatio * (1 - dRatio) / std::pow(dDist, 2);

			for (int nCoord = 0; nCoord < 3; ++nCoord)
			{
				double dProduct = dMagnitude * getSignedDiff(m_vAtoms[nFirst], m_vAtoms[nSecond], m_dBoxSize, nCoord);
				m_vAtoms[nFirst].getForce()[nCoord] += dProduct;
				m_vAtoms[nSecond].getForce()[nCoord] -= dProduct;
			}
		}
	}
}

void Simulation::correctPositions()
{
	for (Atom &cAtom : m_vAtoms)
	{
		std::array<double, 3> adCorrectedPosition = cAtom.getPos();

		for (int nCoord = 0; nCoord < 3; ++nCoord)
		{
			if (adCorrectedPosition[nCoord] < 0.0)
				adCorrectedPosition[nCoord] += m_dBoxSize;
			if (adCorrectedPosition[nCoord] >= m_dBoxSize)
				adCorrectedPosition[nCoord] -= m_dBoxSize;
		}

		cAtom.setPos(adCorrectedPosition);
	}

	// This function doesn't affect periodic distances so no need to update them
}

void Simulation::initialiseDistances()
{
	m_vvdDistances.clear();
	m_vvdDistances.reserve(m_nAtoms);

	for (int nFirst = 0; nFirst < m_nAtoms; ++nFirst)
	{
		m_vvdDistances.emplace_back();
		m_vvdDistances.back().reserve(m_nAtoms);

		for (int nSecond = 0; nSecond < m_nAtoms; ++nSecond)
			if (nFirst >= nSecond)
				m_vvdDistances.back().push_back(0.0f);
			else
				m_vvdDistances.back().push_back(getPeriodicDist(m_vAtoms[nFirst], m_vAtoms[nSecond], m_dBoxSize));
	}
}

void Simulation::updateDistances()
{
	for (int nFirst = 0; nFirst < m_nAtoms; ++nFirst)
	{
		for (int nSecond = nFirst + 1; nSecond < m_nAtoms; ++nSecond)
			m_vvdDistances[nFirst][nSecond] = getPeriodicDist(m_vAtoms[nFirst], m_vAtoms[nSecond], m_dBoxSize);
	}
}
