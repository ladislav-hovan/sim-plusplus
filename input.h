#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <cstring>

#include "params.h"

using std::vector;
using std::array;
using std::string;
using vectorad = vector<array<double, 3> >;

vectorad readValuesFromFile(const string &strFile);

class Input
{
public:
	Input() = default;
	Input(const string& strFilename) { readParameters(strFilename); }
	~Input() = default;

	void readParameters(const string& strFilename);
	void checkValues();

	const InputParams& const getParams() { return m_sParams; }

private:
	InputParams m_sParams{};

	void assignInput(const string& strLine);
};