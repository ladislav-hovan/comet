/*
 * input.h
 *
 *  Created on: 10 Mar 2016
 *      Author: ladislav
 */

#ifndef INPUT_H_
#define INPUT_H_

#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cctype>
#include "small_tools.h"
#include "error_codes.h"
#include "output.h"
#include "input_data.h"

using std::vector;
using std::string;
using vector2d = vector< vector<double> >;

string loadCommandLineInput(int argc, char *argv[]);
InputData loadInput(string strInputFile="parameters.input");
void checkInput(InputData &sInput);
void assignInput(string strLine, InputData &sInput);
vector<double> loadEbetacList(InputData &sInput);
vector<double> loadFesFiles(InputData &sInput);
double getEbetac(string strFilename, InputData &sInput);
vector<double> load1dDataFromFile(string strFilename);
vector2d load2dDataFromFile(string strFilename);

#endif /* INPUT_H_ */
