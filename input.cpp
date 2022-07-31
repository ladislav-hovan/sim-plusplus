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
	if (m_sParams.dBoxSize <= 2 * m_sParams.lj_par.cutoff)
	{
		std::cerr << "The box size needs to be bigger than twice the LJ cutoff distance\n";
		exit(error::invalidInput);
	}
}
