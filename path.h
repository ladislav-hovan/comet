/*
 * path.h
 *
 *  Created on: 17 Mar 2016
 *      Author: ladislav
 */

#ifndef PATH_H_
#define PATH_H_

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <complex>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <limits>
#include <fstream>
#include "input_data.h"
#include "error_codes.h"
#include "output.h"

using std::vector;
using std::string;
using vector2d = vector< vector<double> >;

const double dInf = std::numeric_limits<double>::infinity();

struct PathData
{
	vector<int> vnSnapshots;  // Also includes length via size()
	double dEnergy;
	double dLambda;  // So we don't need to recalculate first differences
	unsigned int nViolations;
};

class Path 
{
public:
	Path();
	Path(vector<double> &vdFes, InputData &sInput, double dBestEnergy);
	Path(vector<double> &vdPathEnergies, vector<double> &vdCoefficients, InputData &sInput,
		 vector2d &vvdNewColvarsPath);
	Path(vector<double> &vdNewGrid, vector<double> &vdFes, vector<double> &vdCoefficients, InputData &sInput,
		 vector2d &vvdNewColvarsPath, bool bAssign=true);
	~Path();

	void optimise();
	double getEnergy();  // Faster way of getting energy, but relies on doing changes via changeAndRecalculate (as optimise does)
	double recomputeEnergy();  // A slower way of getting energy, but doesn't rely on correct use of other functions
	double getBestEnergy();  // Returns energy stored in the BestPath struct
	double getSpectralGap();
	int getBarriers();
	int getViolations();
	double getLambda();
	vector<double> getEigenvalues(int nEigenvalues=-1);
	vector<int> getSnapshots();
	void saveDistances(string strFilename="dist_matrix.dat");

private:
	// Private constructor to be called by other constructors
	Path(vector<double> &vdCoefficients, InputData &sInput, vector2d &vvdNewColvarsPath);

	vector<double> m_vdCoefficients;
	vector2d m_vvdMatrix;
	vector<int> m_vnSnapshots;
	vector<double> m_vdEnergies;
	vector<double> m_vdPeriodicRanges;
	Eigen::MatrixXd m_kMatrix;
	vector<double> m_vdEigenvalues;
	int m_nBarriers = -1;  // Intentionally invalid initial value
	int m_nBarriersSet = -1;  // Same
	double m_dSpectralGap = -1.0;  // Same
	double m_dkT = -1.0;  // Same
	double m_dThreshold = -1.0;  // This technically isn't invalid, but will be overwritten anyway
	double m_dSmoothRatio = -1.0;  // Intentionally invalid
	unsigned int m_nTotalSnapshots = 0;
	unsigned int m_nLengthTarget = 0;
	unsigned int m_nLengthTolerance = 0;
	unsigned int m_nViolations = 0;
	bool m_bStrict = false;
	bool m_bSmoothCount = false;
	bool m_bVerbose = false;

	// Holds best path data
	PathData sBestPath;

	// Extra variables to facilitate optimisation
	vector<double> m_vdFirstDistances;
	vector<double> m_vdSecondDistances;
	double m_dFirstAverage = 0.0;
	double m_dSecondAverage = 0.0;

	// Extra functions for optimisation
	void initialiseOptimisationVars();
	void optimiseSetLength();
	void changeAndRecalculate(unsigned int nPosition, int nSnapshot);
	int checkDistances(unsigned int nPosition);
	int checkAllDistances();
	double calculateDeviation(vector<double> &vdData, double dAverage);
	double calculateDeviation(vector<double> &vdData);

	// Initialisation functions
	void initialiseSnapshots(unsigned int nSnapshots, bool bPrint=true);
	void setPeriodicity(InputData &sInput);
	void loadMatrix(vector<double> &vdCoefficients, vector2d &vvdNewColvarsPath);
	void assignEnergies(vector<double> &vdGrid, vector<double> &vdFes, vector2d &vvdNewColvarsPath);
	void assignEnergies(vector2d &vvdGridList, vector2d &vvdFesList, vector2d &vvdNewColvarsPath);

	// Functions for spectral gap calculation
	void populatekMatrix();
	void calculateSpectralGap();
	void determineSpectralGap();
	int countBarriers(vector<double> &vdFes);
	vector<double> smoothSurface(vector<double> &vdFes);

	// Output functions
	void printPath();
	void logPathEnergies();

	// Friend functions for output
	friend void writeDataFile(string strFilename, vector<double> &vdData);
};

void printLogLine(std::ofstream &Output, InputData &sInput, Path &cPath, vector<double> &vdCoefficients, int nInfCount);

#endif /* PATH_H_ */
