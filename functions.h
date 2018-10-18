/*
 * functions.h
 *
 *  Created on: 18 Mar 2016
 *      Author: ladislav
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <vector>
#include <cmath>
#include <utility>
#include <random>
#include <string>
#include <set>
#include <limits>
#include <fstream>
#include "input.h"
#include "output.h"
#include "path.h"

using std::vector;
using vector2d = vector< vector<double> >;
using pvd = std::pair<vector<double>, vector<double>>;
const double dPi = 3.1415926535;

struct Solution 
{
	vector<double> vdCoefficients;
	vector<double> vdPeriodicCoeffs;
	vector<int> vnSnapshots;
	int nInfinities = -1;
	int nBarriers = -1;
	int nViolations = std::numeric_limits<int>::infinity();
	double dSpectralGap = -1.0;
	double dLambda = 0.0;
};

// Basic math
unsigned int power(unsigned int nBase, unsigned int nExponent);

// Calculation of energy on a grid
vector<bool> determineUpper(unsigned int nPoint, unsigned int nMax);
vector<double> assignSnapshotEnergies(vector2d &vvdFesCompleteData, const vector<int> &vnColvarColumns, 
	int nFreeEnergyColumnComplete, vector2d &vvdColvarPathData);
vector<int> getGridPointLower(vector2d &vvdColvarGrid, vector<double> &vdSpacings, vector<double> &vdCoordinates);
double getGridPointEnergy(vector2d &vvdFesCompleteData, vector<int> &vnLowerPositions, vector<bool> &vbUpper,
	vector2d &vvdColvarGrid, int nEnergyColumn);

// SECTION OF REWEIGHTING FUNCTIONS
// Reweighting functions for 2D surface
void examineNewColvars(vector2d &vvdNewColvars, vector2d &vvdNewGridLimits);
vector2d createNewGrid(vector2d &vvdNewGridLimits, vector<int> &vnGridPoints);
vector<double> calculateNewFes(vector<double> &vdBias, vector2d &vvdNewColvars, vector2d &vvdNewGrid,
    vector<double> &vdEbetacList, vector<int> &vnGridPoints, double dkT, int nFesFilesCount);
std::pair<vector2d, vector<double>> reweightFes(vector<double> &vdBias, vector2d &vvdNewColvars, 
	vector2d &vvdNewColvarsPath, vector<double> &vdEbetacList, InputData &sInput);

// Reweighting with calculation of new variable values
int examineNewColvars(vector< vector<double> > &vvdNewColvars, vector<double> &vdCoefficients,
	std::pair<double, double> &pdNewGridLimits, int nColvars);
vector<double> createNewGrid(std::pair<double, double> pdNewGridLimits, int nGrid);
vector<double> calculateNewFes(vector<double> &vdBias, vector< vector<double> > &vvdNewColvars, 
	vector<double> &vdCoefficients, vector<double> &vdNewGrid, vector<double> &vdEbetacList, vector<double> &vdScalings, 
	InputData &sInput, int nRows);
std::pair<vector<double>, vector<double>> reweightFes(vector<double> &vdBias, vector2d &vvdNewColvars, 
	vector<double> &vdCoefficients, vector<double> &vdEbetacList, vector<double> &vdScalings, InputData &sInput);

// Reweighting without calculation of new variable values but with limiting values (used for Path)
int examineNewColvars(vector<double> &vdNewColvar, std::pair<double, double> &pdNewGridLimits);
vector<double> calculateNewFes(vector<double> &vdBias, vector<double> &vdNewColvar, vector<double> &vdLimits,
	vector<double> &vdNewGrid, vector<double> &vdEbetacList, vector<double> &vdScalings, InputData &sInput, int nRows);
std::pair<vector<double>, vector<double>> reweightFes(vector<double> &vdBias, vector<double> &vdNewColvar, 
	vector<double> &vdLimits, vector<double> &vdEbetacList, vector<double> &vdScalings, InputData &sInput);
// END OF THE SECTION

// Initialisation and manipulation of coefficients vectors
vector<double> initialiseCoefficients(int nColvars);
void normaliseCoefficients(vector<double> &vdCoefficients);
vector<double> perturbCoefficients(vector<double> &vdCoefficients, std::mt19937 &Mersenne);
vector<double> perturbPeriodicCoefficients(vector<double> &vdCoefficients, std::mt19937 &Mersenne, 
	vector<double> &vdPeriodicRanges);
vector<double> randomiseCoefficients(int nNumber, std::mt19937 &Mersenne);
vector<double> randomisePeriodicCoefficients(int nNumber, std::mt19937 &Mersenne, vector<double> &vdPeriodicRanges);

// Replacement of infinities in the free energy surface
int correctInfinities(vector<double> &vdNewFes, bool bAltInfinity);
void replaceInfinities(vector<double> &vdNewFes);
void interpolateInfinities(vector<double> &vdNewFes);

// Various auxiliary functions
vector2d convertPeriodic(vector2d &vvdColvars, vector<double> &vdPeriodicCoeffs, InputData &sInput);
void decreasekT(InputData &sInput);
void setMinToZero(vector<double> &vdEnergies);

// Path collective variable calculation
void calculatePath(vector2d &vvdNewColvarsPath, vector2d &vvdNewColvars, vector<double> vdCoefficients, 
	vector<int> vnSnapshots, double dLambda, InputData &sInput, pvd &pvdPathValues);

// Master functions and logging
void optimiseAndLog(Path &cPath, std::ofstream &LogOutput, Solution &sCurrent, InputData &sInput);
Solution stopTrying(Solution &sCurrent, double dLambda, Path &cPath1, vector<int> &vnSnapshots, 
	vector<double> &vdCoefficients, int nCountInfs, double dBestEnergy, std::ofstream &LogOutput);
Solution tryCoefficients(vector<double> &vdCoefficients, vector<double> &vdPeriodicCoeffs, vector<double> &vdSummedBias,
	vector2d &vvdNewColvars, vector<double> &vdEbetacList, vector<double> &vdScalings, InputData &sInput, 
	std::ofstream &LogOutput);
Solution tryCoefficientsPath(vector<double> &vdCoefficients, vector<double> &vdPeriodicCoeffs, vector<double> &vdSummedBias,
	vector2d &vvdNewColvars, vector2d &vvdNewColvarsPath, vector<double> &vdEbetacList, vector<double> &vdScalings, 
	InputData &sInput, std::ofstream &LogOutput);

#endif /* FUNCTIONS_H_ */
