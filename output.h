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
#include <experimental/filesystem>
#include "input_data.h"

using std::vector;
namespace fs = std::experimental::filesystem;

void writeDataFile(std::string strFilename, vector<double> &vdData);
void writeDataFile(std::string strFilename, vector< vector<double> > &vvdData);
void backupFile(std::string strFilename);
void printNumbersToFile(std::ofstream &Output, std::string strDecoder, ...);
void createColvarFile(vector< vector<double> > &vvdNewColvarsPath, vector<int> &vnSnapshots, std::string strFilename="COLVAR_PATH");

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
