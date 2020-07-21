/*
 * auxiliary.cpp
 *
 *  Created on: 16 Mar 2016
 *      Author: ladislav
 */

#include "auxiliary.h"

vector<double> sumColumns(vector2d &vvdData, vector<int> &vnColumns) 
{
	vector<double> vdSummedData;

	for (unsigned int nRow = 0; nRow < vvdData.size(); ++nRow)
	{
		double dSum = 0.0;

		for (unsigned int nColumn = 0; nColumn < vnColumns.size(); ++nColumn)
		{
			dSum += vvdData.at(nRow).at(vnColumns.at(nColumn));
		}

		vdSummedData.push_back(dSum);
	}

	return vdSummedData;
}

vector<double> sumColumns(vector2d &vvdData) 
{
	vector<double> vdSummedData;

	for (unsigned int nRow = 0; nRow < vvdData.size(); ++nRow)
	{
		double dSum = 0.0;

		for (unsigned int nColumn = 0; nColumn < vvdData.at(nRow).size(); ++nColumn)
		{
			dSum += vvdData.at(nRow).at(nColumn);
		}

		vdSummedData.push_back(dSum);
	}

	return vdSummedData;
}

vector2d selectColumns(vector2d &vvdData, vector<int> &vnColumns, InputData &sInput, bool bPath) 
{
	vector2d vvdSelectedData;

	static bool bScaled = false;
	static bool bScaledPath = false;

	if (bPath && !bScaledPath)
	{
		for (auto &nCol : sInput.vnPeriodicColumnsPath)
			nCol = static_cast<int>(std::find(vnColumns.begin(), vnColumns.end(), nCol) - vnColumns.begin());
		bScaledPath = true;
	}
	else if (!bScaled)
	{
		for (auto &nCol : sInput.vnPeriodicColumns)
			nCol = static_cast<int>(std::find(vnColumns.begin(), vnColumns.end(), nCol) - vnColumns.begin());
		bScaled = true;
	}

	for (unsigned int nRow = 0; nRow < vvdData.size(); ++nRow)
	{
		vector<double> vdSelection;

		unsigned int nNext = 0;
		for (unsigned nColumn = 0; nColumn < vvdData.at(nRow).size(); ++nColumn)
		{
			if (nColumn == vnColumns.at(nNext))
			{
				vdSelection.push_back(vvdData.at(nRow).at(nColumn));
				if (nNext < vnColumns.size() - 1)
					++nNext;
				else
					break;
			}
		}

		vvdSelectedData.push_back(vdSelection);
	}

	return vvdSelectedData;
}

vector2d getLimits(vector2d &vvdDataLimits)
{
	vector2d vvdLimits;

	for (unsigned int nRow = 0; nRow < vvdDataLimits.size(); ++nRow)
	{
		if (vvdLimits.empty())
		{
			vector<double> vdLimits;

			for (unsigned int nColumn = 0; nColumn < vvdDataLimits.at(nRow).size(); ++nColumn)
				vdLimits.push_back(vvdDataLimits.at(nRow).at(nColumn));

			// Push back twice for lower and upper limit
			vvdLimits.push_back(vdLimits);
			vvdLimits.push_back(vdLimits);
		}
		else
		{
			for (unsigned int nColumn = 0; nColumn < vvdDataLimits.at(nRow).size(); ++nColumn)
			{
				double dValue = vvdDataLimits.at(nRow).at(nColumn);

				if (dValue < vvdLimits.at(0).at(nColumn))
					vvdLimits.at(0).at(nColumn) = dValue;
				if (dValue > vvdLimits.at(1).at(nColumn))
					vvdLimits.at(1).at(nColumn) = dValue;
			}
		}
	}

	return vvdLimits;
}

vector2d rescaleDataRange(vector2d &vvdDataLimits, vector2d &vvdDataToScale, InputData &sInput, bool bPath) 
{
	vector2d vvdScaledData;

	if (vvdDataLimits.at(0).size() != vvdDataToScale.at(0).size())
	{
		std::cerr << "The number of columns in the files for rescaling does not match" << std::endl;
		exit (SIZE_NOT_MATCHED);
	}

	vector2d vvdLimits = getLimits(vvdDataLimits);

	vector<double> vdFactors;
	for (unsigned int nColumn = 0; nColumn < vvdLimits.at(0).size(); ++nColumn)
		vdFactors.push_back(vvdLimits.at(1).at(nColumn) - vvdLimits.at(0).at(nColumn));

	// Static so only rescaled once, even if this function is invoked multiple times
	static bool bPeriodicScaled = false;
	if (!bPeriodicScaled)
	{
		for (unsigned int nCount = 0; nCount < sInput.vdPeriodicRanges.size(); ++nCount)
		{
			if (bPath)
				sInput.vdPeriodicRanges.at(nCount) /= vdFactors.at(sInput.vnPeriodicColumnsPath.at(nCount));
			else
				sInput.vdPeriodicRanges.at(nCount) /= vdFactors.at(sInput.vnPeriodicColumns.at(nCount));
		}
		bPeriodicScaled = true;
	}

	for (unsigned int nRow = 0; nRow < vvdDataToScale.size(); ++nRow)
	{
		vector<double> vdRescaledRow;

		for (unsigned int nColumn = 0; nColumn < vvdDataToScale.at(nRow).size(); ++nColumn)
		{
			double dRescaled = (vvdDataToScale.at(nRow).at(nColumn) - vvdLimits.at(0).at(nColumn)) / vdFactors.at(nColumn);
			vdRescaledRow.push_back(dRescaled);
		}

		vvdScaledData.push_back(vdRescaledRow);
	}

	return vvdScaledData;
}

vector2d createTestingList(int nCoefficients) 
{
	vector2d vvdTestingList;

	for (int nRow = 0; nRow < nCoefficients; ++nRow)
	{
		vector<double> vdRow;
		for (int nColumn = 0; nColumn < nCoefficients; ++nColumn)
		{
			if (nRow == nColumn)
				vdRow.push_back(1.0);
			else
				vdRow.push_back(0.0);
		}
		vvdTestingList.push_back(vdRow);
	}

	return vvdTestingList;
}
