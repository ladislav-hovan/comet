/*
 * inputdata.h
 *
 *  Created on: 12 Apr 2016
 *      Author: ladislav
 */

#ifndef INPUTDATA_H_
#define INPUTDATA_H_

#include <vector>
#include <string>

using std::vector;
using std::string;

struct InputData 
{
	// The files to take input from
	string strMetadColvarFile;
	string strNewColvarFile;
	string strNewColvarFilePath;
	string strFesPrefix;
	string strFesDataFile;
	string strFesComplete;

	// Further specifications of the files
	int nScalingColumn = -1;
	int nFreeEnergyColumn = -1;
	int nFesFilesCount = -1;
	int nFreeEnergyColumnComplete = -1;
	vector<int> vnBiasColumns;
	vector<int> vnNewColumns;
	vector<int> vnNewColumnsPath;
	vector<int> vnColvarColumns;
	vector<int> vnColvarColumnsPath;
	string strColvarFilePath;

	// Periodicity specifications
	vector<int> vnPeriodicColumns;
	vector<int> vnPeriodicColumnsPath;
	vector<double> vdPeriodicRanges;

	// Additional specifications
	double dkT = -1.0;
	double dGamma = -1.0;  // Biasfactor, 0 for non-tempered metadynamics
	int nGrid = -1;  // Number of bins to reweight to
	int nLengthTarget = -1;  // Number of snapshots to use
	double dAnnealingkT = -1.0;  // Initial kT for simulated annealing
	int nAnnealingSteps = -1;

	// Optional arguments
	string strLogFile = "log_spectral_gaps.dat";
	bool bEuclidean = false;
	bool bAltInfinity = false;
	bool bLogEnergy = false;
	int nLogEigenvalues = 0;
	bool bRescale = true;
	int nBarriers = -1;
	int nTimeColumn = 0;
	string strCoeffList;
	bool bTryAll = false;
	bool bRandom = false;
	int nLengthTolerance = 0;
	bool bStrictOrder = false;
	bool bSPath = false;
	bool bPath2D = false;
	bool bLimitIncorrect = false;
	double dInfLimit = 0.0;
	double dCooling = -1.0;
	double dThreshold = 1.0;
	double dZLimit = -1.0;
	bool bConvertPeriodic = false;
	string strPeriodicCoeffList;
	bool bSmoothCount = false;
	bool bNoPath = false;
	int nSeed = -1;
};

#endif /* INPUTDATA_H_ */
