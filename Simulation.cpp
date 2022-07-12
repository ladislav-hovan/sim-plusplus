#include "Simulation.h"

Simulation::Simulation(InputParams& sInput): m_dTimeStep{sInput.dTimeStep}, m_dBoxSize{sInput.dBoxSize}, m_nSeed{sInput.nSeed}, m_nMaxSteps{sInput.nSteps}, m_nAtoms{sInput.nAtoms},
											 m_dMass{sInput.dMass}, m_dTemp{sInput.dTemp}, m_LJPar{sInput.lj_par}
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
		std::array<double, 3> adNewPosition = cAtom.getPos();
		std::array<double, 3> adVelocity = cAtom.getVelocity();
		std::array<double, 3> adOldForce = cAtom.getOldForce();

		for (int nCoord = 0; nCoord < 3; ++nCoord)
		{
			adNewPosition[nCoord] += adVelocity[nCoord] * m_dTimeStep;
			adNewPosition[nCoord] += adOldForce[nCoord] * std::pow(m_dTimeStep, 2) / (2 * cAtom.getMass());
		}

		cAtom.setPos(adNewPosition);

		// TODO: Add some check to warn against big changes in position
	}
}

void Simulation::updateVelocities()
{
	for (Atom &cAtom : m_vAtoms)
	{
		std::array<double, 3> adNewVelocity = cAtom.getVelocity();
		std::array<double, 3> adOldForce = cAtom.getOldForce();
		std::array<double, 3> adForce = cAtom.getForce();

		for (int nCoord = 0; nCoord < 3; ++nCoord)
			adNewVelocity[nCoord] += (adOldForce[nCoord] + adForce[nCoord]) * m_dTimeStep / (2 * cAtom.getMass());

		cAtom.setVelocity(adNewVelocity);

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
			double dDist = getPeriodicDist(m_vAtoms[nFirst], m_vAtoms[nSecond], m_dBoxSize);
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
}