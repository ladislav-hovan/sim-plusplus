#include "input.h"

bool isEmpty(string &strInput)
{
	int nCount = 0;

	while (strInput[nCount])
	{
		if (!std::isspace(strInput[nCount])) 
			return false;
		++nCount;
	}

	return true;
}

string trimWhitespace(const string &strInput)
{
	size_t nFirst = strInput.find_first_not_of(" \t");
	size_t nSecond = strInput.find_last_not_of(" \t");

	string strOutput = strInput.substr(nFirst, nSecond - nFirst + 1);

	return strOutput;
}

vector<double> separateStringDoubles(const string &strSequence)
{
	vector<double> vdSequence;
	string strMember;
	int nCount = 0;

	while (strSequence[nCount])
	{
		if (!strncmp(&strSequence[nCount], ",", 1))
		{
			double dMember = std::stod(strMember);
			vdSequence.push_back(dMember);
			strMember.clear();
		}
		else
			strMember += strSequence[nCount];

		++nCount;
	}

	double dMember = std::stod(strMember);
	vdSequence.push_back(dMember);
	strMember.clear();

	return vdSequence;
}

vectorad readValuesFromFile(const string& strFile)
{
	vectorad vadValues;

	std::ifstream InputStream(strFile);

	if (!InputStream)
	{
		std::cerr << "Couldn't open file " << strFile << std::endl;
		exit(error::fileNotFound);
	}

	while (InputStream)
	{
		string strLine;
		getline(InputStream, strLine);

		// Ignore empty lines
		if (isEmpty(strLine))
			continue;

		array<double, 3> adTemp{0, 0, 0};

		std::stringstream LineStream(strLine);

		for (int nPos = 0; nPos < 3; ++nPos)
			LineStream >> adTemp[nPos];

		vadValues.push_back(adTemp);
	}

	return vadValues;
}

void Input::assignInput(const string& strLine)
{
	string strFirstInput, strSecondInput;

	size_t nPos = strLine.find_first_of('=');
	if (nPos == string::npos)
	{
		std::cerr << "The following line doesn't contain any assignment:" << std::endl;
		std::cerr << strLine << std::endl;
		exit(error::noAssignment);
	}
	else
	{
		strFirstInput = trimWhitespace(strLine.substr(0, nPos));
		strSecondInput = trimWhitespace(strLine.substr(nPos + 1));
	}

	// Assigning strings, simple
	if (!strFirstInput.compare("EnergyFile"))			m_sParams.strEnergyFile = strSecondInput;
	if (!strFirstInput.compare("PositionFile"))			m_sParams.strPositionFile = strSecondInput;
	if (!strFirstInput.compare("PosInputFile"))			m_sParams.strPosInputFile = strSecondInput;
	if (!strFirstInput.compare("VelInputFile"))			m_sParams.strVelInputFile = strSecondInput;
	if (!strFirstInput.compare("BoxType"))				m_sParams.strBoxType = strSecondInput;

	// Assigning integers
	if (!strFirstInput.compare("Seed"))					m_sParams.nSeed = std::stoi(strSecondInput);
	if (!strFirstInput.compare("MaxSteps"))				m_sParams.nSteps = std::stoi(strSecondInput);
	if (!strFirstInput.compare("NumAtoms"))				m_sParams.nAtoms = std::stoi(strSecondInput);
	if (!strFirstInput.compare("InitialStep"))			m_sParams.nInitialStep = std::stoi(strSecondInput);

	// Assigning doubles
	if (!strFirstInput.compare("LJ_Epsilon"))			m_sParams.lj_par.epsilon = std::stod(strSecondInput);
	if (!strFirstInput.compare("LJ_Cutoff"))			m_sParams.lj_par.cutoff = std::stod(strSecondInput);
	if (!strFirstInput.compare("LJ_Rmin"))				m_sParams.lj_par.r_m = std::stod(strSecondInput);
	if (!strFirstInput.compare("Temperature"))			m_sParams.dTemp = std::stod(strSecondInput);
	if (!strFirstInput.compare("BoxSize"))				m_sParams.dBoxSize = std::stod(strSecondInput);
	if (!strFirstInput.compare("Mass"))					m_sParams.dMass = std::stod(strSecondInput);
	if (!strFirstInput.compare("TimeStep"))				m_sParams.dTimeStep = std::stod(strSecondInput);
	if (!strFirstInput.compare("InitialTime"))			m_sParams.dInitialTime = std::stod(strSecondInput);

	// Assigning vectors of doubles
	if (!strFirstInput.compare("BoxVectors"))			m_sParams.vdBoxVectors = separateStringDoubles(strSecondInput);
}

void Input::readParameters(const string& strFilename)
{
	std::ifstream InputStream(strFilename);

	if (!InputStream)
	{
		std::cerr << "Couldn't open file " << strFilename << std::endl;
		exit(error::fileNotFound);
	}

	while (InputStream)
	{
		string strLine;
		getline(InputStream, strLine);

		// Ignore comment or empty lines
		if (!std::strncmp(strLine.c_str(), "#", 1) || isEmpty(strLine))
			continue;

		// Remove everything after #
		size_t nPos = strLine.find_first_of('#');
		if (nPos != string::npos)
			strLine.erase(nPos);

		assignInput(strLine);
	}
}

void Input::checkValues()
{
	// Check the parameters for consistency
	if ((m_sParams.vdBoxVectors.size() > 0) && (m_sParams.dBoxSize > 0.0))
	{
		std::cerr << "Both box vectors and box size were specified in the input file\n";
		exit(error::invalidInput);
	}

	if ((m_sParams.vdBoxVectors.size() > 0) && (!m_sParams.strBoxType.compare("cubic")))
	{
		std::cerr << "Box vectors instead of box size are being specified for a cubic box\n";
		exit(error::invalidInput);
	}

	if (m_sParams.vdBoxVectors.size() == 9)
	{
		// The first three correspond to v1(x), v2(y) and v3(z)

		// TODO: Implement vector rotation to make sure v1(y) = v1(z) = v2(z) = 0
	}
	else if (m_sParams.vdBoxVectors.size() == 3)
	{
		// Fill in the other members with zero
		while (m_sParams.vdBoxVectors.size() < 9)
			m_sParams.vdBoxVectors.push_back(0.0f);
	}
	else if (m_sParams.vdBoxVectors.size() == 0)
	{
		if (m_sParams.dBoxSize < 0)
		{
			// None of them are defined, use the default value
			if (!m_sParams.strBoxType.compare("cubic"))
			{
				// Box is cubic
				m_sParams.dBoxSize = 10.0f;  // nm
			}
			else
			{
				std::cerr << "A non-cubic box was selected but no size parameters were provided\n";
				exit(error::invalidInput);
			}
		}
		else
		{
			if (m_sParams.strBoxType.compare("cubic"))
			{
				std::cerr << "A non-cubic box was selected but the box vectors were not provided\n";
				exit(error::invalidInput);
			}
		}
	}
	else
	{
		std::cerr << "The length of the provided box vectors needs to be 9 or 3\n";
		exit(error::invalidInput);
	}

	if (!m_sParams.strBoxType.compare("cubic"))
	{
		if (m_sParams.dBoxSize <= 2 * m_sParams.lj_par.cutoff)
		{
			std::cerr << "The box size needs to be bigger than twice the LJ cutoff distance\n";
			exit(error::invalidInput);
		}
	}
}
