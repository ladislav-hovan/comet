#pragma once

#include <vector>
#include <limits>

using std::vector;

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