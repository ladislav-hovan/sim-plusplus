#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <cstring>

#include "input.h"
#include "params.h"

using std::vector;
using std::array;
using std::string;
using vectorad = vector<array<double, 3> >;

InputParams readParameters(const string &strFilename);
vectorad readValuesFromFile(const string &strFile);