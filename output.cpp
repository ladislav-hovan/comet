/*
 * output.cpp
 *
 *  Created on: 16 Mar 2016
 *      Author: ladislav
 */

#include "output.h"

void writeDataFile(std::string strFilename, vector<double> &vdData)
{
	std::ofstream OutputStream(strFilename.c_str());

	for (unsigned int nRow = 0; nRow < vdData.size(); ++nRow)
	{
		OutputStream << vdData.at(nRow) << "\n";
	}
}

void writeDataFile(std::string strFilename, vector< vector<double> > &vvdData) 
{
	std::ofstream OutputStream(strFilename.c_str());

	for (unsigned int nRow = 0; nRow < vvdData.size(); ++nRow)
	{
		for (unsigned int nColumn = 0; nColumn < vvdData.at(nRow).size(); ++nColumn)
		{
			OutputStream << vvdData.at(nRow).at(nColumn) << " ";
		}

		OutputStream << "\n";
	}
}

void backupFile(std::string strFilename) 
{
	fs::path Path = strFilename;

	if (fs::exists(Path))
	{
		int nBackup = 0;
		fs::path PathNew;
		do
		{
			std::stringstream Stream;
			Stream << "bck." << nBackup << "." << strFilename;
			PathNew = Stream.str();
			++nBackup;
		}
		while (fs::exists(PathNew));

		std::cout << "Backing up file " << Path << " to " << PathNew << std::endl;
		fs::rename(Path, PathNew);
	}
}

void createColvarFile(vector< vector<double> > &vvdNewColvarsPath, vector<int> &vnSnapshots, std::string strFilename) 
{
	std::cout << "Creating a new COLVAR file called: " << strFilename << std::endl;
	backupFile(strFilename);

	std::ofstream OutputStream(strFilename.c_str());
	for (unsigned int nSnapshot = 0; nSnapshot < vnSnapshots.size(); ++nSnapshot)
	{
		OutputStream << nSnapshot + 1 << " ";

		int nRow = vnSnapshots[nSnapshot];
		for (unsigned int nColumn = 0; nColumn < vvdNewColvarsPath[nRow].size(); ++nColumn)
			OutputStream << vvdNewColvarsPath[nRow][nColumn] << " ";
		OutputStream << "\n";
	}
}
