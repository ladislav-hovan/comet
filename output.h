/*
 * output.h
 *
 *  Created on: 16 Mar 2016
 *      Author: ladislav
 */

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdarg>
//#include <experimental/filesystem>
#include <filesystem>
#include "input_data.h"
#include "solution.h"

using std::vector;
using std::string;
using vector2d = vector< vector<double> >;

//namespace fs = std::experimental::filesystem;
namespace fs = std::filesystem;

void writeDataFile(string strFilename, vector<double> &vdData);
void writeDataFile(string strFilename, vector2d &vvdData);
void backupFile(string strFilename);
void createColvarFile(vector2d &vvdNewColvarsPath, vector<int> &vnSnapshots, int nPrecision, 
	string strFilename="COLVAR_PATH");
void createPlumedFile(Solution &sBest, vector2d &vvdLimits, InputData &sInput, string strFilename="plumed_path.dat");

// Template functions
template <typename T>
void printVector(vector<T> &vtData, bool bTerminate=false) 
{
	for (unsigned int nCount = 0; nCount < vtData.size(); ++nCount)
	{
		std::cout << vtData.at(nCount);
		if (nCount < vtData.size() - 1)
			std::cout << ", ";
	}

	if (bTerminate)
		std::cout << std::endl;
}

template <typename T>
void printNumberToFile(std::ofstream &Output, T tValue, bool bTerminate=false) 
{
	Output << tValue << " ";
	if (bTerminate)
		Output << "\n";
}

template <typename T>
void printNumbersToFile(std::ofstream &Output, vector<T> vtValues, bool bTerminate=false)
{
	for (auto tValue: vtValues)
		Output << tValue << " ";
	if (bTerminate)
		Output << "\n";
}

#endif /* OUTPUT_H_ */
