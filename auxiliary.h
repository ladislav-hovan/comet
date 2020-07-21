/*
 * auxiliary.h
 *
 *  Created on: 16 Mar 2016
 *      Author: ladislav
 */

#ifndef AUXILIARY_H_
#define AUXILIARY_H_

#include <vector>
#include <iostream>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "error_codes.h"
#include "input_data.h"

using std::vector;
using vector2d = vector< vector<double> >;

vector<double> sumColumns(vector2d &vvdData, vector<int> &vnColumns);
vector<double> sumColumns(vector2d &vvdData);
vector2d selectColumns(vector2d &vvdData, vector<int> &vnColumns, InputData &sInput, bool bPath);
vector2d getLimits(vector2d &vvdDataLimits);
vector2d rescaleDataRange(vector2d &vvdDataLimits, vector2d &vvdDataToScale, InputData &sInput, bool bPath);
vector2d createTestingList(int nCoefficients);

// Template functions
template <typename T>
vector<T> reduceRows(vector<T> &vtData, int nTarget) 
{
	vector<T> vtReduced;

	int nCurrent = vtData.size();
	double dRatio = static_cast<double>(nCurrent) / nTarget;
	double dTargetCount = 0.0;

	for (int nCount = 0; nCount < nCurrent; ++nCount)
	{
		if (nCount >= dTargetCount)
		{
			vtReduced.push_back(vtData.at(nCount));
			dTargetCount += dRatio;
		}
	}

	return vtReduced;
}

template <typename T>
void cutExtraTimes(vector< vector<T> > &vvtData, double dEndTime, int nTimeColumn) 
{
	double dSmallestDist = std::numeric_limits<double>::infinity();
	int nIndex = -1;

	for (unsigned int nRow = 0; nRow < vvtData.size(); ++nRow)
	{
		double dDiff = static_cast<double>(vvtData[nRow][nTimeColumn]) - dEndTime;
		double dDist = fabs(dDiff);

		if (dDist < dSmallestDist)
		{
			dSmallestDist = dDist;
			nIndex = nRow;
		}
	}
	vvtData.erase(vvtData.begin() + nIndex + 1, vvtData.end());
}

template <typename T1, typename T2>
void alignEndTime(vector< vector<T1> > &vvtData1, vector< vector<T2> > &vvtData2, int nTimeColumn) 
{
	double dEndTime1 = static_cast<double>(vvtData1[vvtData1.size() - 1][nTimeColumn]);
	double dEndTime2 = static_cast<double>(vvtData2[vvtData2.size() - 1][nTimeColumn]);

	if (dEndTime1 < dEndTime2)
		cutExtraTimes(vvtData2, dEndTime1, nTimeColumn);
	else if (dEndTime1 > dEndTime2)
		cutExtraTimes(vvtData1, dEndTime2, nTimeColumn);
	// If they are equal, no alignment needed
}

template <typename T1, typename T2>
void alignRows(vector<T1> &vtData1, vector<T2> &vtData2) 
{
	int nSize1 = vtData1.size();
	int nSize2 = vtData2.size();

	if (nSize1 == nSize2)
		return;
	else if (nSize1 > nSize2)
		vtData1 = reduceRows(vtData1, nSize2);
	else  // nSize2 > nSize1
		vtData2 = reduceRows(vtData2, nSize1);
}

template <typename T1, typename T2>
void alignTimes(vector< vector<T1> > &vvtDataLong, vector< vector<T2> > &vvtDataShort, int nTimeColumn) 
{
	auto Target = vvtDataShort.begin();
	auto tTargetTime = (*Target)[nTimeColumn];
	vector< vector<T1> > vvtShortened;
	bool bDone = false;
	T1 tPreviousTime = -1;  // Intentionally invalid

	for (unsigned int nCount = 0; nCount < vvtDataLong.size(); ++nCount)
	{
		if (nCount == 0)
		{
			tPreviousTime = vvtDataLong[0][nTimeColumn];

			while (tPreviousTime >= tTargetTime)
			{
				vvtShortened.push_back(vvtDataLong[0]);
				++Target;

				if (Target == vvtDataShort.end())
				{
					bDone = true;
					break;
				}

				tTargetTime = (*Target)[nTimeColumn];
			}
		}

		auto tCurrentTime = vvtDataLong[nCount][nTimeColumn];

		while (!bDone && tTargetTime < tCurrentTime)
		{
			if (tCurrentTime > tTargetTime && tPreviousTime <= tTargetTime)
			{
				if ((tCurrentTime - tTargetTime) < (tTargetTime - tPreviousTime))
					vvtShortened.push_back(vvtDataLong[nCount]);
				else
					vvtShortened.push_back(vvtDataLong[nCount - 1]);

				++Target;

				if (Target == vvtDataShort.end())
				{
					bDone = true;
					break;
				}

				tTargetTime = (*Target)[nTimeColumn];
			}
		}
		if (bDone) 
			break;

		tPreviousTime = tCurrentTime;
	}

	vvtDataLong = vvtShortened;
}

template <typename T1, typename T2>
void alignRows(vector< vector<T1> > &vvtData1, vector< vector<T2> > &vvtData2, int nTimeColumn) 
{
	if (nTimeColumn < 0)
		// Length alignment only
		alignRows(vvtData1, vvtData2);
	else
	{
		alignEndTime(vvtData1, vvtData2, nTimeColumn);

		int nSize1 = vvtData1.size();
		int nSize2 = vvtData2.size();

		if (nSize1 > nSize2)
			alignTimes(vvtData1, vvtData2, nTimeColumn);
		else if (nSize1 < nSize2)
			alignTimes(vvtData2, vvtData1, nTimeColumn);
		// Else we're aligned already
	}
}

#endif /* AUXILIARY_H_ */
