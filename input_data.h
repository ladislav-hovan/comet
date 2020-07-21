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
	int nFreeEnergyColumn = -1;
	int nFesFilesCount = -1;
	int nFreeEnergyColumnComplete = -1;
	vector<int> vnBiasColumns;
	vector<int> vnRBiasColumns;
	vector<int> vnNewColumns;
	vector<int> vnNewColumnsPath;
	vector<int> vnColvarColumns;
	vector<int> vnColvarColumnsPath;
	string strColvarFilePath;
	vector<double> vdLowerLimits;
	vector<double> vdUpperLimits;

	// Periodicity specifications
	vector<int> vnPeriodicColumns;
	vector<int> vnPeriodicColumnsPath;
	vector<double> vdPeriodicRanges;

	// Additional specifications
	double dkT = -1.0;
	double dGamma = -1.0;  // Biasfactor, 0 for non-tempered metadynamics
	int nGrid = -1;  // Number of bins to reweight to
	int nLengthTarget = -1;  // Number of snapshots to use
	int nLengthTolerance = 0;
	double dAnnealingkT = -1.0;  // Initial kT for simulated annealing
	int nAnnealingSteps = -1;

	// Optional arguments
	// Strings
	string strLogFile = "log_spectral_gaps.dat";
	string strCoeffList;
	string strPeriodicCoeffList;
	// Integers
	int nLogEigenvalues = 0;
	int nBarriers = -2;
	int nTimeColumn = 0;
	int nSeed = -1;
	int nMinBarriers = 0;
	int nPrecision = 5;
	// Doubles
	double dInfLimit = 0.0;
	double dCooling = -1.0;
	double dThreshold = 1.0;
	double dZLimit = -1.0;
	double dPerturbationLimit = 0.3;
	double dSmoothRatio = 0.03;
	// Booleans
	bool bAltInfinity = true;
	bool bLogEnergy = true;
	bool bSmoothCount = true;
	bool bRescale = false;
	bool bDescale = false;
	bool bTryAll = false;
	bool bRandom = false;
	bool bStrictOrder = false;
	bool bSPath = false;
	bool bPath2D = false;
	bool bLimitIncorrect = false;
	bool bConvertPeriodic = false;
	bool bNoPath = false;
	bool bVerbose = false;
};

#endif /* INPUTDATA_H_ */
