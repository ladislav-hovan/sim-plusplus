#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "params.h"

using std::string;

InputParams readParameters(const string &strFilename);