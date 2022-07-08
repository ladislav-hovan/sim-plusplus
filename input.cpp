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

string trimWhitespace(string &strInput)
{
	size_t nFirst = strInput.find_first_not_of(" \t");
	size_t nSecond = strInput.find_last_not_of(" \t");

	string strOutput = strInput.substr(nFirst, nSecond - nFirst + 1);

	return strOutput;
}

void assignInput(string &strLine, InputParams &sParams)
{
	string strFirstInput, strSecondInput;

	size_t nPos = strLine.find_first_of('=');
	if (nPos == string::npos)
	{
		std::cerr << "The following line doesn't contain any assignment:" << std::endl;
		std::cerr << strLine << std::endl;
		exit(3);
	}
	else
	{
		strFirstInput = trimWhitespace(strLine.substr(0, nPos));
		strSecondInput = trimWhitespace(strLine.substr(nPos + 1));
	}

	// Assigning strings, simple
	if (!strFirstInput.compare("EnergyFile"))			sParams.strEnergyFile = strSecondInput;
	if (!strFirstInput.compare("PositionFile"))			sParams.strPositionFile = strSecondInput;

	// Assigning integers
	if (!strFirstInput.compare("Seed"))					sParams.nSeed = std::stoi(strSecondInput);
	if (!strFirstInput.compare("MaxSteps"))				sParams.nSteps = std::stoi(strSecondInput);
	if (!strFirstInput.compare("NumAtoms"))				sParams.nAtoms = std::stoi(strSecondInput);

	// Assigning doubles
	if (!strFirstInput.compare("LJ_Epsilon"))			sParams.lj_par.epsilon = std::stod(strSecondInput);
	if (!strFirstInput.compare("LJ_Cutoff"))			sParams.lj_par.cutoff = std::stod(strSecondInput);
	if (!strFirstInput.compare("LJ_Rmin"))				sParams.lj_par.r_m = std::stod(strSecondInput);
	if (!strFirstInput.compare("Temperature"))			sParams.dTemp = std::stod(strSecondInput);
	if (!strFirstInput.compare("BoxSize"))				sParams.dBoxSize = std::stod(strSecondInput);
	if (!strFirstInput.compare("Mass"))					sParams.dMass = std::stod(strSecondInput);
	if (!strFirstInput.compare("TimeStep"))				sParams.dTimeStep = std::stod(strSecondInput);
}

InputParams readParameters(const string &strFilename)
{
	InputParams sParams;
	
	std::ifstream InputStream(strFilename);
	
	if (!InputStream)
	{
		std::cerr << "Couldn't open file " << strFilename << std::endl;
		exit(2);
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

		assignInput(strLine, sParams);
	}

	return sParams;
}
