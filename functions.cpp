/*
 * functions.cpp
 *
 *  Created on: 18 Mar 2016
 *      Author: ladislav
 */

#include "functions.h"

unsigned int power(unsigned int nBase, unsigned int nExponent) 
{
	int nResult = 1;

	for (unsigned int nCount = 0; nCount < nExponent; ++nCount)
		nResult *= nBase;

	return nResult;
}

vector<bool> determineUpper(unsigned int nPoint, unsigned int nMax) 
{
	vector<bool> vbUpper;

	for (unsigned int nPos = 0; nPos < nMax; ++nPos)
	{
		vbUpper.push_back(nPoint % 2);
		nPoint /= 2;  // Integer division intentional
	}

	return vbUpper;
}

double getGridPointEnergy(vector2d &vvdFesCompleteData, vector<int> &vnLowerPositions, vector<bool> &vbUpper,
	vector2d &vvdColvarGrid, int nEnergyColumn) 
{
	vector<int> vnPositions = vnLowerPositions;

	for (unsigned int nCount = 0; nCount < vbUpper.size(); ++nCount)
	{
		if (vbUpper[nCount])
			++vnPositions[nCount];
	}

	int nRow = vnPositions[0];
	int nToMultiply = 1;

	for (unsigned int nCount = 1; nCount < vvdColvarGrid.size(); ++nCount)
	{
		nToMultiply *= vvdColvarGrid[nCount - 1].size();
		nRow += vnPositions[nCount] * nToMultiply;
	}

	return vvdFesCompleteData[nRow][nEnergyColumn];
}

vector<int> getGridPointLower(vector2d &vvdColvarGrid, vector<double> &vdSpacings, vector<double> &vdCoordinates) 
{
	vector<int> vnLowerElements;

	for (unsigned int nCoord = 0; nCoord < vdCoordinates.size(); ++nCoord)
	{
		double dPosition = (vdCoordinates[nCoord] - vvdColvarGrid[nCoord][0]) / vdSpacings[nCoord];
		vnLowerElements.push_back(static_cast<int>(floor(dPosition)));
	}

	return vnLowerElements;
}

vector<double> assignSnapshotEnergies(vector2d &vvdFesCompleteData, const vector<int> &vnColvarColumns, int nFreeEnergyColumnComplete,
									  vector2d &vvdColvarPathData) 
{
	// Determine grid values for every CV from the FES
	// Use set, so they don't repeat
	unsigned int nMax = vnColvarColumns.size();
	vector<std::set<double>> vsColvarGrid(nMax);

	for (auto vdFesRow: vvdFesCompleteData)
	{
		for (unsigned int nColumn = 0; nColumn < nMax; ++nColumn)
		{
			double dColvarValue = vdFesRow[vnColvarColumns[nColumn]];
			vsColvarGrid[nColumn].insert(dColvarValue);
		}
	}

	// Convert to vector for fast access, log spacings too
	vector2d vvdColvarGrid;
	vector<double> vdSpacings;

	for (auto sColvarRange: vsColvarGrid)
	{
		vector<double> vdColvarValues;

		for (auto dValue: sColvarRange)
			vdColvarValues.push_back(dValue);

		vvdColvarGrid.push_back(vdColvarValues);
		vdSpacings.push_back(vdColvarValues[1] - vdColvarValues[0]);
	}

	// For every snapshot, assign the most likely value by interpolating
	vector<double> vdSnapshotEnergies;

	for (auto vdColvarDataRow: vvdColvarPathData)
	{
		vector<int> vnLowerPositions = getGridPointLower(vvdColvarGrid, vdSpacings, vdColvarDataRow);
		double dResult = 0.0;

		for (unsigned int nPoint = 0; nPoint < power(2, vnLowerPositions.size()); ++nPoint)
		{
			vector<bool> vbUpper = determineUpper(nPoint, nMax);
			double dContribution = getGridPointEnergy(vvdFesCompleteData, vnLowerPositions, vbUpper, vvdColvarGrid,
				nFreeEnergyColumnComplete);

			for (unsigned int nCount = 0; nCount < vbUpper.size(); ++nCount)
			{
				if (vbUpper[nCount])
					dContribution *= (vvdColvarGrid[nCount][vnLowerPositions[nCount] + 1] - vdColvarDataRow[nCount]);
				else
					dContribution *= (vdColvarDataRow[nCount] - vvdColvarGrid[nCount][vnLowerPositions[nCount]]);

				dContribution /= vdSpacings[nCount];
			}

			dResult += dContribution;
		}

		vdSnapshotEnergies.push_back(dResult);
	}

	return vdSnapshotEnergies;
}

void setMinToZero(vector<double> &vdEnergies) 
{
	double dMin = dInf;

	for (const double &dValue: vdEnergies)
		if (dValue < dMin)
			dMin = dValue;

	for (double &dValue: vdEnergies)
		dValue -= dMin;
}

// New functions to facilitate 2D reweighting
void examineNewColvars(vector2d &vvdNewColvars, vector2d &vvdNewGridLimits)
{
	for (auto &row: vvdNewColvars)
	{
		bool bSkip = false;

		if (vvdNewGridLimits.empty())
		{
			for (double &val: row)
				if (std::isinf(val))
					bSkip = true;

			if (bSkip) 
				continue;

			vvdNewGridLimits.push_back(row);
			vvdNewGridLimits.push_back(row);
		}
		else
		{
			for (double &val: row)
				if (std::isinf(val))
					bSkip = true;

			if (bSkip) 
				continue;

            for (unsigned int nCol = 0; nCol < row.size(); ++nCol)
            {
                if (row[nCol] < vvdNewGridLimits[0][nCol])
                    vvdNewGridLimits[0][nCol] = row[nCol];
                if (row[nCol] > vvdNewGridLimits[1][nCol])
                    vvdNewGridLimits[1][nCol] = row[nCol];
            }
		}
	}

	// Add extra 0.1% to range, to avoid values at the very edge
	for (unsigned int nCol = 0; nCol < vvdNewGridLimits[0].size(); ++nCol)
	{
		double dRange = vvdNewGridLimits[1][nCol] - vvdNewGridLimits[0][nCol];

		vvdNewGridLimits[0][nCol] -= 0.0005 * dRange;
		vvdNewGridLimits[1][nCol] += 0.0005 * dRange;
	}
}

vector2d createNewGrid(vector2d &vvdNewGridLimits, vector<int> &vnGridPoints)
{
    vector2d vvdNewGrid;
    int nSize = 1;

    for (int nPoints: vnGridPoints)
        nSize *= nPoints;

    vvdNewGrid.reserve(nSize);
    vector<double> vdDeltas;
    int nColvars = vvdNewGridLimits[0].size();

    for (int nColvar = 0; nColvar < nColvars; ++nColvar)
        vdDeltas.push_back((vvdNewGridLimits[1][nColvar] - vvdNewGridLimits[0][nColvar])
                           / (vnGridPoints[nColvar] - 1));

    vector<int> vnCounter(vnGridPoints.size(), 0);

    do
    {
        vector<double> vdColvarRow;

        for (int nColvar = 0; nColvar < nColvars; ++nColvar)
            vdColvarRow.push_back(vvdNewGridLimits[0][nColvar] + vnCounter[nColvar] * vdDeltas[nColvar]);

        vvdNewGrid.push_back(vdColvarRow);
        ++vnCounter.front();

        for (int nColvar = 0; nColvar < nColvars; ++nColvar)
        {
            if (vnCounter[nColvar] == vnGridPoints[nColvar])
            {
                if (nColvar < nColvars - 1)
                {
                    vnCounter[nColvar] = 0;
                    vnCounter[nColvar + 1] += 1;
                }
            }
        }
    }
    while (vnCounter.back() != vnGridPoints.back());

    return vvdNewGrid;
}

vector<double> calculateNewFes(vector<double> &vdBias, vector2d &vvdNewColvars, vector2d &vvdNewGrid,
                               vector<double> &vdEbetacList, vector<int> &vnGridPoints, double dkT, int nFesFilesCount)
{
    int nColvars = vvdNewGrid[0].size();
    int nRows = vvdNewColvars.size();
    int nGrid = vvdNewGrid.size();

    // Initialise the FES to zero
    vector<double> vdFes(nGrid, 0.0);

    double dDenom = 0.0;
    double dRatio = double(nFesFilesCount) / nRows;
    vector<double> vdCutoffs;

    for (int nColvar = 0; nColvar < nColvars; ++nColvar)
        vdCutoffs.push_back((vvdNewGrid.back()[nColvar] - vvdNewGrid[0][nColvar])
                            / (2 * (vnGridPoints[nColvar] - 1)));

    printVector(vdCutoffs, true);

    for (int nRow = 0; nRow < nRows; ++nRow)
	{
        double dBias = vdBias[nRow];
        int nLocation = 0, nColvar = nColvars - 1;

        for (int nPos = 0; nPos < nGrid; ++nPos)
        {
            bool bDone = false;

            while (fabs(vvdNewGrid[nPos][nColvar] - vvdNewColvars[nRow][nColvar]) < vdCutoffs[nColvar])
            {
                if (nColvar > 0)
                    --nColvar;
                else
                {
                    nLocation = nPos;
                    bDone = true;
                    break;
                }
            }

            if (bDone)
                break;
        }

        // Assumes continuous time progression of work estimates and colvar rows
        int nIndex = int(ceil(float(nRow + 1) * dRatio)) - 1;

        double dExpBias = exp(dBias / dkT) / vdEbetacList[nIndex];
        vdFes[nLocation] += dExpBias;
        dDenom += dExpBias;
    }

    // Calculate FES values, find the minimum
    double dMin = dInf;

    for (int nPos = 0; nPos < nGrid; ++nPos)
    {
        vdFes[nPos] = -dkT * log(vdFes[nPos] / dDenom);

        if (vdFes[nPos] < dMin)
            dMin = vdFes[nPos];
    }

    // Set minimum to zero
    for (double &dValue: vdFes)
        dValue -= dMin;

    return vdFes;
}

std::pair<vector2d, vector<double>> reweightFes(vector<double> &vdBias, vector2d &vvdNewColvars, 
	vector2d &vvdNewColvarsPath, vector<double> &vdEbetacList, InputData &sInput)
{
	vector2d vvdNewGridLimits;
	examineNewColvars(vvdNewColvars, vvdNewGridLimits);
	examineNewColvars(vvdNewColvarsPath, vvdNewGridLimits);

	printVector(vvdNewGridLimits[0], true);
	printVector(vvdNewGridLimits[1], true);

	vector<int> vnGridPoints(2, sInput.nGrid);
	vector2d vvdNewGrid = createNewGrid(vvdNewGridLimits, vnGridPoints);

	vector<double> vdFes = calculateNewFes(vdBias, vvdNewColvars, vvdNewGrid, vdEbetacList, vnGridPoints, sInput.dkT, sInput.nFesFilesCount);

	return std::pair<vector2d, vector<double>> (vvdNewGrid, vdFes);
}
// End of 2D functions

vector<double> calculateNewFes(vector<double> &vdBias, vector2d &vvdNewColvars, vector<double> &vdCoefficients,
	vector<double> &vdNewGrid, vector<double> &vdEbetacList, vector<double> &vdScalings, InputData &sInput, int nRows) 
{
	int nGrid = sInput.nGrid;
	int nColvars = sInput.vnNewColumns.size();
	double dkT = sInput.dkT;

    // Initialise the FES to zero
    vector<double> vdFes(nGrid, 0.0);

    double dDenom = 0.0;
    double dRatio = double(sInput.nFesFilesCount) / nRows;
    double dCutoff = (vdNewGrid[1] - vdNewGrid[0]) / 2;

    for (int nRow = 0; nRow < nRows; ++nRow)
    {
    	double dBias = vdBias[nRow];
        vector<double> *p_vdColvarRow = &(vvdNewColvars[nRow]);
        double dTotalValue = 0.0;

        for (int nColumn = 0; nColumn < nColvars; ++nColumn)
        	dTotalValue += (*p_vdColvarRow)[nColumn] * vdCoefficients[nColumn];

		// Using the work estimates
		if (vdScalings.empty())
		{
			int nLocation = 0;

			for (int nPos = 0; nPos < nGrid; ++nPos)
			{
				double dDiff = fabs(vdNewGrid[nPos] - dTotalValue);

				if (dDiff < dCutoff)
				{
					nLocation = nPos;
					break;
				}
			}

			int nIndex = int(ceil(float(nRow + 1) * dRatio)) - 1;

			double dExpBias = exp(dBias / dkT) / vdEbetacList[nIndex];
			vdFes[nLocation] += dExpBias;
			dDenom += dExpBias;
		}
		// Using the weights to create a weighted histogram
		else
		{
			double dOverflow = 0.5 * (1 + std::erf((dTotalValue - vdNewGrid[0] / (2 * dCutoff)))
				+ (1 + std::erfc((dTotalValue - vdNewGrid[nGrid] / (2 * dCutoff)))));
			double dScaling = vdScalings[nRow] / (1.0 - dOverflow);

			for (int nPos = 0; nPos < nGrid; ++nPos)
			{
				double dDiff = vdNewGrid[nPos] - dTotalValue;
				double dKernel = (1 / sqrt(2 * dPi * dCutoff)) * exp(-pow(dDiff, 2) / (2 * pow(dCutoff, 2)));
				vdFes[nPos] += dScaling * dKernel;
			}

			dDenom += dScaling;
		}
    }

    // Calculate FES values, find the minimum
    double dMin = dInf;

    for (int nPos = 0; nPos < nGrid; ++nPos)
    {
        vdFes[nPos] = - dkT * log(vdFes[nPos] / dDenom);

        if (vdFes[nPos] < dMin)
			dMin = vdFes[nPos];
    }

    // Set minimum to zero
    for (double &dValue: vdFes)
    	dValue -= dMin;

    return vdFes;
}

vector<double> createNewGrid(std::pair<double, double> pdNewGridLimits, int nGrid) 
{
	vector<double> vdNewGrid;

	double dDelta = (pdNewGridLimits.second - pdNewGridLimits.first) / (nGrid - 1);

	for (int nCount = 0; nCount < nGrid; ++nCount)
		vdNewGrid.push_back(pdNewGridLimits.first + dDelta * nCount);

	return vdNewGrid;
}

int examineNewColvars(vector< vector<double> > &vvdNewColvars, vector<double> &vdCoefficients,
						std::pair<double, double> &pdNewGridLimits, int nColvars) 
{
	int nRows = vvdNewColvars.size();

	for (int nCount = 0; nCount < nRows; ++nCount)
	{
		double dTotalValue = 0.0;

		for (int nColvar = 0; nColvar < nColvars; ++nColvar)
			dTotalValue += vdCoefficients.at(nColvar) * vvdNewColvars.at(nCount).at(nColvar);

		if (dTotalValue < pdNewGridLimits.first) 
			pdNewGridLimits.first = dTotalValue;
		if (dTotalValue > pdNewGridLimits.second) 
			pdNewGridLimits.second = dTotalValue;
	}

	return nRows;
}

pvd reweightFes(vector<double> &vdBias, vector2d &vvdNewColvars, vector<double> &vdCoefficients, vector<double> &vdEbetacList, 
	vector<double> &vdScalings, InputData &sInput) 
{
	std::pair<double, double> pdNewGridLimits = {dInf, -dInf};
	int nRows = examineNewColvars(vvdNewColvars, vdCoefficients, pdNewGridLimits, sInput.vnNewColumns.size());

	vector<double> vdNewGrid = createNewGrid(pdNewGridLimits, sInput.nGrid);

	vector<double> vdFes = calculateNewFes(vdBias, vvdNewColvars, vdCoefficients, vdNewGrid, vdEbetacList, vdScalings, sInput, nRows);

	return pvd (vdNewGrid, vdFes);
}

int examineNewColvars(vector<double> &vdNewColvar, std::pair<double, double> &pdNewGridLimits) 
{
	int nRows = vdNewColvar.size();

	for (double dNewValue: vdNewColvar)
	{
		if (dNewValue < pdNewGridLimits.first) pdNewGridLimits.first = dNewValue;
		if (dNewValue > pdNewGridLimits.second) pdNewGridLimits.second = dNewValue;
	}

	return nRows;
}

vector<double> calculateNewFes(vector<double> &vdBias, vector<double> &vdNewColvar, vector<double> &vdLimits,
	vector<double> &vdNewGrid, vector<double> &vdEbetacList, vector<double> &vdScalings, InputData &sInput, int nRows) 
{
	int nGrid = sInput.nGrid;
	int nColvars = sInput.vnNewColumns.size();
	double dkT = sInput.dkT;

    // Initialise the FES to zero
    vector<double> vdFes(nGrid, 0.0);

    double dDenom = 0.0;
    double dRatio = double(sInput.nFesFilesCount) / nRows;
    double dCutoff = (vdNewGrid[1] - vdNewGrid[0]) / 2;

    // Testing applying z limits only to endpoints
    double dMinS = vdNewGrid[0] + 0.5;
    double dMaxS = vdNewGrid[vdNewGrid.size() - 1] - 0.5;

    for (int nRow = 0; nRow < nRows; ++nRow)
    {
        // Depends on the order of evaluation
        if (vdLimits.empty() || ((vdNewColvar[nRow] <  dMinS || vdNewColvar[nRow] > dMaxS) && (vdLimits[nRow] < sInput.dZLimit)) ||
            ((dMinS <= vdNewColvar[nRow]) && (vdNewColvar[nRow] <= dMaxS)))
        {
            double dBias = vdBias[nRow];
            double dTotalValue = vdNewColvar[nRow];

			// Using the work estimates
			if (vdScalings.empty())
			{
				int nLocation = 0;

				for (int nPos = 0; nPos < nGrid; ++nPos)
				{
					double dDiff = fabs(vdNewGrid[nPos] - dTotalValue);

					if (dDiff < dCutoff)
					{
						nLocation = nPos;
						break;
					}
				}

				int nIndex = int(ceil(float(nRow + 1) * dRatio)) - 1;

				double dExpBias = exp(dBias / dkT) / vdEbetacList[nIndex];
				vdFes[nLocation] += dExpBias;
				dDenom += dExpBias;
			}
			// Using the weights to create a weighted histogram
			else
			{
				double dScaling = vdScalings[nRow];

				for (int nPos = 0; nPos < nGrid; ++nPos)
				{
					double dDiff = vdNewGrid[nPos] - dTotalValue;
					double dKernel = (1 / sqrt(2 * dPi * dCutoff)) * exp(-pow(dDiff, 2) / (2 * pow(dCutoff, 2)));
					vdFes[nPos] += dScaling * dKernel;
				}

				dDenom += dScaling;
			}
        }
    }

    // Calculate FES values, find the minimum
    double dMin = dInf;

    for (int nPos = 0; nPos < nGrid; ++nPos)
    {
		vdFes[nPos] = - dkT * log(vdFes[nPos] / dDenom);

        if (vdFes[nPos] < dMin) 
			dMin = vdFes[nPos];
    }

    // Set minimum to zero
    for (double &dValue: vdFes)
    	dValue -= dMin;

    return vdFes;
}

std::pair<vector<double>, vector<double>> reweightFes(vector<double> &vdBias, vector<double> &vdNewColvar, vector<double> &vdLimits,
	vector<double> &vdEbetacList, vector<double> &vdScalings, InputData &sInput) 
{
	std::pair<double, double> pdNewGridLimits = {dInf, -dInf};
	int nRows = examineNewColvars(vdNewColvar, pdNewGridLimits);

	vector<double> vdNewGrid = createNewGrid(pdNewGridLimits, sInput.nGrid);

	vector<double> vdFes = calculateNewFes(vdBias, vdNewColvar, vdLimits, vdNewGrid, vdEbetacList, vdScalings, sInput, nRows);

	return pvd (vdNewGrid, vdFes);
}

vector<double> initialiseCoefficients(int nColvars) 
{
	vector<double> vdCoefficients(nColvars, 0.0);

	for (double &dCoefficient: vdCoefficients)
		dCoefficient = sqrt(1.0 / nColvars);

	return vdCoefficients;
}

void normaliseCoefficients(vector<double> &vdCoefficients) 
{
	double dLength = 0.0;

	for (const double &dValue: vdCoefficients)
		dLength += pow(dValue, 2);

	dLength = sqrt(dLength);

	for (double &dValue: vdCoefficients)
		dValue /= dLength;
}

vector<double> perturbCoefficients(vector<double> &vdCoefficients, std::mt19937 &Mersenne) 
{
	vector<double> vdNewCoefficients = vdCoefficients;
	std::normal_distribution<double> Distribution(0.0, 0.3);
	double dLength = 0.0;

	for (double &dValue: vdNewCoefficients)
	{
		dValue *= (1.0 + Distribution(Mersenne));

		while (dValue < 0.0)
			dValue += 1.0;
		while (dValue > 1.0)
			dValue -= 1.0;

		dLength += pow(dValue, 2);
	}

	dLength = sqrt(dLength);

	// Normalise
	for (double &dValue: vdNewCoefficients)
		dValue /= dLength;

	return vdNewCoefficients;
}

vector<double> perturbPeriodicCoefficients(vector<double> &vdCoefficients, std::mt19937 &Mersenne, vector<double> &vdPeriodicRanges) 
{
	vector<double> vdNewCoefficients = vdCoefficients;
	std::normal_distribution<double> Distribution(0.0, 0.3);

	for (unsigned int nCount = 0; nCount < vdNewCoefficients.size(); ++nCount)
	{
		// Testing all the same
		if (nCount == 0)
		{
			double dLimit = vdPeriodicRanges[nCount] / 2;
			double &dValue = vdNewCoefficients[nCount];
			dValue += Distribution(Mersenne);

			while (dValue < 0)
				dValue += dLimit;
			while (dValue > dLimit)
				dValue -= dLimit;
		}
		else
			vdNewCoefficients[nCount] = vdNewCoefficients[nCount - 1];
	}

	return vdNewCoefficients;
}

vector<double> randomiseCoefficients(int nNumber, std::mt19937 &Mersenne) 
{
	vector<double> vdCoefficients(nNumber, 0.0);
	std::uniform_real_distribution<double> Distribution(0.0, 1.0);
	double dLength = 0.0;

	for (double &dValue: vdCoefficients)
	{
		dValue = Distribution(Mersenne);
		dLength += pow(dValue, 2);
	}

	dLength = sqrt(dLength);

	// Normalise
	for (double &dValue: vdCoefficients)
		dValue /= dLength;

	return vdCoefficients;
}

vector<double> randomisePeriodicCoefficients(int nNumber, std::mt19937 &Mersenne, vector<double> &vdPeriodicRanges) 
{
	vector<double> vdCoefficients(nNumber, 0.0);

	for (int nCount = 0; nCount < nNumber; ++nCount)
	{
		double dLimit = vdPeriodicRanges[nCount] / 2;
		std::uniform_real_distribution<double> Distribution(-dLimit, dLimit);
		vdCoefficients[nCount] = Distribution(Mersenne);
	}

	return vdCoefficients;
}

int correctInfinities(vector<double> &vdNewFes, bool bAltInfinity) 
{
	int nInfCount = 0;

	for (const auto &vdValue: vdNewFes)
	{
		if (std::isinf(vdValue))
			++nInfCount;
	}

	if (nInfCount > 0)
	{
		if (nInfCount == 1)
			std::cout << "Detected 1 infinity in the reweighted FES" << std::endl;
		else
			std::cout << "Detected " << nInfCount << " infinities in the reweighted FES" << std::endl;

		if (bAltInfinity)
			replaceInfinities(vdNewFes);
		else
			interpolateInfinities(vdNewFes);
	}
	else
		std::cout << "No infinities detected in the reweighted FES" << std::endl;

	return nInfCount;
}

void replaceInfinities(vector<double> &vdNewFes) 
{
	double dFesMax = -dInf;

	for (const double &dValue: vdNewFes)
	{
		if (!std::isinf(dValue))
		{
			if (dValue > dFesMax)
				dFesMax = dValue;
		}
	}

	for (double &dValue: vdNewFes)
	{
		if (std::isinf(dValue))
			dValue = dFesMax;
	}
}

void interpolateInfinities(vector<double> &vdNewFes)
{
	int nPrevious = -1;
	int nInfinities = 0;

	for (unsigned int nCount = 0; nCount < vdNewFes.size(); ++nCount)
	{
		if (std::isinf(vdNewFes[nCount]))
		{
			++nInfinities;

			if (nCount == vdNewFes.size() - 1)  // Last value infinite, shouldn't happen but to be sure
			{
				double dGradient = vdNewFes[nPrevious] - vdNewFes[nPrevious - 1];
				for (unsigned int nPos = nPrevious + 1, nFactor = 1; nPos <= nCount; ++nPos, ++nFactor)
					vdNewFes[nPos] = vdNewFes[nPrevious] + dGradient * nFactor;
			}
		}
		else
		{
			if (nInfinities > 0)
			{
				if (nPrevious == -1)
				{
					if (std::isinf(vdNewFes[nCount + 1]))  // First finite value, followed by infinity again
					{
						for (unsigned int nPos = 0; nPos < nCount; ++nPos)
							vdNewFes[nPos] = vdNewFes[nCount];
					}
					else  // First finite value, followed by another one
					{
						double dGradient = vdNewFes[nCount + 1] - vdNewFes[nCount];
						for (unsigned int nPos = 0; nPos < nCount; ++nPos)
							vdNewFes[nPos] = vdNewFes[nCount] + dGradient * (nInfinities - nPos);
					}
				}
				else
				{
					double dGradient = (vdNewFes[nCount] - vdNewFes[nPrevious]) / (nCount - nPrevious);
					for (unsigned int nPos = nPrevious + 1, nFactor = 1; nPos < nCount; ++nPos, ++nFactor)
						vdNewFes[nPos] = vdNewFes[nPrevious] + dGradient * nFactor;
				}
				nInfinities = 0;
			}
			nPrevious = nCount;
		}
	}
}

void decreasekT(InputData &sInput) 
{
	static double dDeltakT = sInput.dAnnealingkT / sInput.nAnnealingSteps;

	if (sInput.dCooling < 0.0)
		sInput.dAnnealingkT -= dDeltakT;
	else
		sInput.dAnnealingkT *= sInput.dCooling;
}

vector2d convertPeriodic(vector2d &vvdColvars, vector<double> &vdPeriodicCoeffs, InputData &sInput) 
{
	vector2d vvdConverted;
	static vector<double> vdPeriodics;

	if (vdPeriodics.empty())
	{
		for (unsigned int nColumn = 0; nColumn < sInput.vnNewColumnsPath.size(); ++nColumn)
		{
			int nPeriodic = -1;
			for (unsigned int nCount = 0; nCount < sInput.vnPeriodicColumnsPath.size(); ++nCount)
				if (nColumn == sInput.vnPeriodicColumnsPath[nCount])
					nPeriodic = nCount;

			if (nPeriodic < 0)
				vdPeriodics.push_back(-1.0);
			else
				vdPeriodics.push_back(sInput.vdPeriodicRanges[nPeriodic]);
		}
	}

	for (unsigned int nRow = 0; nRow < vvdColvars.size(); ++nRow)
	{
		vector<double> vdRow;
		int nCount = 0;

		for (unsigned int nColumn = 0; nColumn < vvdColvars[nRow].size(); ++nColumn)
		{
			if (vdPeriodics[nColumn] >= 0.0)
			{
				vdRow.push_back(0.5 + cos(vvdColvars[nRow][nColumn] - vdPeriodicCoeffs[nCount]));
				++nCount;
			}
			else
				vdRow.push_back(vvdColvars[nRow][nColumn]);
		}

		vvdConverted.push_back(vdRow);
	}

	return vvdConverted;
}

void calculatePath(vector2d &vvdNewColvarsPath, vector2d &vvdNewColvars, vector<double> vdCoefficients, vector<int> vnSnapshots,
				   double dLambda, InputData &sInput, pvd &pvdPathValues) 
{
	std::cout << "Calculating the path collective variables" << std::endl;

	static vector<double> vdPeriodicRanges;

	if (vdPeriodicRanges.empty())
	{
		if (sInput.bConvertPeriodic)
		{
			for (unsigned int nColumn = 0; nColumn < sInput.vnNewColumnsPath.size(); ++nColumn)
				vdPeriodicRanges.push_back(-1.0);
		}
		else
		{
			for (unsigned int nColumn = 0; nColumn < sInput.vnNewColumnsPath.size(); ++nColumn)
			{
				int nPeriodic = -1;

				for (unsigned int nCount = 0; nCount < sInput.vnPeriodicColumnsPath.size(); ++nCount)
					if (nColumn == sInput.vnPeriodicColumnsPath[nCount])
						nPeriodic = nCount;

				if (nPeriodic < 0)
					vdPeriodicRanges.push_back(-1.0);
				else
					vdPeriodicRanges.push_back(sInput.vdPeriodicRanges[nPeriodic]);
			}
		}
	}

	for (unsigned int nCount = 0; nCount < vvdNewColvars.size(); ++nCount)
	{
		double dSValue = 0.0;
		double dPartition = 0.0;

		for (unsigned int nPoint = 0; nPoint < vnSnapshots.size(); ++nPoint)
		{
			int nRow = vnSnapshots[nPoint];
			double dSum = 0.0;

			for (unsigned int nColvar = 0; nColvar < sInput.vnNewColumnsPath.size(); ++nColvar)
			{
				double dDiff = vvdNewColvars[nCount][nColvar] - vvdNewColvarsPath[nRow][nColvar];

				if (vdPeriodicRanges[nColvar] > 0)
				{
					if (dDiff > (vdPeriodicRanges[nColvar] / 2))
						dDiff -= vdPeriodicRanges[nColvar];
					if (dDiff < (- vdPeriodicRanges[nColvar] / 2))
						dDiff += vdPeriodicRanges[nColvar];
				}

				dSum += pow(dDiff * vdCoefficients[nColvar], 2);
			}

			double dExp = exp(- dLambda * dSum);
			dSValue += (nPoint + 1) * dExp;
			dPartition += dExp;
		}

		dSValue /= dPartition;
		pvdPathValues.first.push_back(dSValue);

		if ((sInput.dZLimit > 0.0) || sInput.bPath2D || sInput.bLimitIncorrect)
			pvdPathValues.second.push_back(-log(dPartition) / dLambda);  // Definition of z variable
	}
}

Solution stopTrying(Solution &sCurrent, double dLambda, Path &cPath1, vector<int> &vnSnapshots, vector<double> &vdCoefficients,
					int nCountInfs, double dBestEnergy, std::ofstream &LogOutput)
{
	sCurrent.dLambda = dLambda;
	sCurrent.nViolations = cPath1.getViolations();
	sCurrent.vnSnapshots = vnSnapshots;

	printNumbersToFile(LogOutput, vdCoefficients);
	printNumberToFile(LogOutput, -1);
	printNumberToFile(LogOutput, -1.0);
	printNumberToFile(LogOutput, nCountInfs);
	printNumberToFile(LogOutput, dBestEnergy);
	LogOutput << "\n";

	return sCurrent;
}

Solution tryCoefficientsPath(vector<double> &vdCoefficients, vector<double> &vdPeriodicCoeffs, vector<double> &vdSummedBias,
	vector2d &vvdNewColvars, vector2d &vvdNewColvarsPath, vector<double> &vdEbetacList, vector<double> &vdScalings, 
	InputData &sInput, std::ofstream &LogOutput)
{
	Solution sCurrent;
	sCurrent.vdCoefficients = vdCoefficients;

	// The reweighting version
	vector2d vvdConverted, vvdConvertedPath;

	// Before anything, convert periodic variables to non-periodic ones
	if (sInput.bConvertPeriodic)
	{
		vvdConverted = convertPeriodic(vvdNewColvars, vdPeriodicCoeffs, sInput);
		vvdConvertedPath = convertPeriodic(vvdNewColvarsPath, vdPeriodicCoeffs, sInput);
		sCurrent.vdPeriodicCoeffs = vdPeriodicCoeffs;

		// writeDataFile("converted_colvar.dat", vvdConvertedPath);
	}

	vector2d &rvvdConverted = sInput.bConvertPeriodic ? vvdConverted : vvdNewColvars;
	vector2d &rvvdConvertedPath = sInput.bConvertPeriodic ? vvdConvertedPath : vvdNewColvarsPath;

	// First, need optimally spaced snapshots but energy is irrelevant
	vector<double> vdEmpty;
	Path cPath1(vdEmpty, vdEmpty, vdCoefficients, sInput, rvvdConvertedPath, false);
	cPath1.optimise();
	vector<int> vnSnapshots = cPath1.getSnapshots();
	double dLambda = cPath1.getLambda();
	double dBestEnergy = cPath1.getBestEnergy();
	cPath1.saveDistances();

	// Second, calculate path variables for every point in new colvars and also path if 2D
	vector<double> vdFirst, vdSecond;
	vdFirst.reserve(vvdNewColvars.size());

	if (sInput.dZLimit > 0.0)
		vdSecond.reserve(vvdNewColvars.size());

	pvd pvdPathValues (vdFirst, vdSecond);
	calculatePath(rvvdConvertedPath, rvvdConverted, vdCoefficients, vnSnapshots, dLambda, sInput, pvdPathValues);
	vector<double> &vdSValues = pvdPathValues.first;
	vector<double> &vdZValues = pvdPathValues.second;
	pvd pvdSnapshotValues;
	vector2d vvdSnapshotPathValues;

	if (sInput.bPath2D)
	{
		calculatePath(rvvdConvertedPath, rvvdConvertedPath, vdCoefficients, vnSnapshots, dLambda, sInput, pvdSnapshotValues);
		printVector(pvdSnapshotValues.first, true);
		printVector(pvdSnapshotValues.second, true);

		for (unsigned int nCount = 0; nCount < pvdSnapshotValues.first.size(); ++nCount)
		{
			vvdSnapshotPathValues.push_back({pvdSnapshotValues.first[nCount], pvdSnapshotValues.second[nCount]});

			if (std::isinf(pvdSnapshotValues.second[nCount]) || std::isnan(pvdSnapshotValues.second[nCount]) ||
				std::isinf(pvdSnapshotValues.first[nCount]) || std::isnan(pvdSnapshotValues.first[nCount]))
				return stopTrying(sCurrent, dLambda, cPath1, vnSnapshots, vdCoefficients, -1, dBestEnergy, LogOutput);
		}

		writeDataFile("snapshots_path_values.dat", vvdSnapshotPathValues);
	}

	// Third, reweight onto the s variable (possibly subject to z limits) or also z variable
	vector<double> vdNewFes, vdSnapshotEnergies;
	vector2d vvdPathValues;

	if (sInput.bLimitIncorrect || sInput.bPath2D)
	{
		// Counting infinities in PCVs and creating vector2d
		int nCountInfs = 0;

		for (unsigned int nCount = 0; nCount < vdSValues.size(); ++nCount)
		{
			if (std::isinf(vdZValues[nCount]) || std::isnan(vdZValues[nCount]) ||
				std::isinf(vdSValues[nCount]) || std::isnan(vdSValues[nCount]))
				++nCountInfs;
			vvdPathValues.push_back({vdSValues[nCount], vdZValues[nCount]});
		}

		std::cout << nCountInfs << " " << vdZValues.size() << std::endl;

		static vector<double> vdInfRatios;
		double dInfRatio = static_cast<double>(nCountInfs) / vdZValues.size();
		vdInfRatios.push_back(dInfRatio);

		// Not sure about this
		// writeDataFile("inf_ratios.dat", vdInfRatios);

		if (dInfRatio > sInput.dInfLimit)
			return stopTrying(sCurrent, dLambda, cPath1, vnSnapshots, vdCoefficients, nCountInfs, dBestEnergy, LogOutput);
	}

	if (sInput.bPath2D)
	{
		std::pair<vector2d, vector<double>> pvdNewFesData = reweightFes(vdSummedBias, vvdPathValues, vvdSnapshotPathValues,
			vdEbetacList, sInput);

		// Outputting the FES and also creating vector for assigning energies to snapshots
		std::ofstream FesTest("generated_fes_2d.dat");
		double dPastZ = pvdNewFesData.first[0][1];
		vector2d vvdNewFesData(pvdNewFesData.first.size(), vector<double>(3));

		for (unsigned int nRow = 0; nRow < pvdNewFesData.first.size(); ++nRow)
		{
			if (pvdNewFesData.first[nRow][1] != dPastZ)
			{
				FesTest << "\n";
				dPastZ = pvdNewFesData.first[nRow][1];
			}

			vvdNewFesData[nRow][0] = pvdNewFesData.first[nRow][0];
			vvdNewFesData[nRow][1] = pvdNewFesData.first[nRow][1];
			vvdNewFesData[nRow][2] = pvdNewFesData.second[nRow];

			printNumbersToFile(FesTest, vvdNewFesData[nRow], true);
		}

		FesTest.close();

		// Assign energies to path snapshots
		vdSnapshotEnergies = assignSnapshotEnergies(vvdNewFesData, {0, 1}, 2, vvdSnapshotPathValues);
		printVector(vdSnapshotEnergies, true);
	}
	else
	{
		pvd pvdNewFesData = reweightFes(vdSummedBias, vdSValues, vdZValues, vdEbetacList, vdScalings, sInput);
		vdNewFes = pvdNewFesData.second;
		writeDataFile("new_fes.dat", vdNewFes);
	}

	if (sInput.bPath2D)
	{
		sCurrent.nInfinities = correctInfinities(vdSnapshotEnergies, sInput.bAltInfinity);

		if (sCurrent.nInfinities > 0)
			return stopTrying(sCurrent, dLambda, cPath1, vnSnapshots, vdCoefficients, sCurrent.nInfinities, dBestEnergy, 
				LogOutput);
	}
	else
		sCurrent.nInfinities = correctInfinities(vdNewFes, sInput.bAltInfinity);

	// Fourth, construct a path corresponding to the FES in s variable and get spectral gap
	Path cPath2;

	if (sInput.bPath2D)
	{
		// Just testing
		cPath2 = Path(vdSnapshotEnergies, vdCoefficients, sInput, vvdNewColvarsPath);
		cPath2.optimise();
	}
	else
		cPath2 = Path(vdNewFes, sInput, dBestEnergy);

	std::cout << "Barriers: " << cPath2.getBarriers() << ", Spectral gap: " << cPath2.getSpectralGap() << std::endl;

	sCurrent.nBarriers = cPath2.getBarriers();
	sCurrent.dSpectralGap = cPath2.getSpectralGap();
	sCurrent.dLambda = dLambda;
	sCurrent.nViolations = cPath1.getViolations();
	sCurrent.vnSnapshots = vnSnapshots;

	printLogLine(LogOutput, sInput, cPath2, sCurrent.vdCoefficients, sCurrent.nInfinities);

	return sCurrent;
}

Solution tryCoefficients(vector<double> &vdCoefficients, vector<double> &vdPeriodicCoeffs, vector<double> &vdSummedBias,
	vector2d &vvdNewColvars, vector<double> &vdEbetacList, vector<double> &vdScalings, InputData &sInput, std::ofstream &LogOutput)
{
	Solution sCurrent;
	sCurrent.vdCoefficients = vdCoefficients;

	// The reweighting version
	vector2d vvdConverted;

	// Before anything, convert periodic variables to non-periodic ones
	if (sInput.bConvertPeriodic)
	{
		vvdConverted = convertPeriodic(vvdNewColvars, vdPeriodicCoeffs, sInput);
		sCurrent.vdPeriodicCoeffs = vdPeriodicCoeffs;

		// writeDataFile("converted_colvar.dat", vvdConvertedPath);
	}

	vector2d &rvvdConverted = sInput.bConvertPeriodic ? vvdConverted : vvdNewColvars;

	pvd pvdNewFesData = reweightFes(vdSummedBias, rvvdConverted, vdCoefficients, vdEbetacList, vdScalings, sInput);
	vector<double> vdNewGrid = pvdNewFesData.first;
	vector<double> vdNewFes = pvdNewFesData.second;
	writeDataFile("new_fes.dat", vdNewFes);

	sCurrent.nInfinities = correctInfinities(vdNewFes, sInput.bAltInfinity);

	Path cPath(vdNewFes, sInput, -1.0);

	std::cout << "Barriers: " << cPath.getBarriers() << ", Spectral gap: " << cPath.getSpectralGap() << std::endl;
	sCurrent.nBarriers = cPath.getBarriers();
	sCurrent.dSpectralGap = cPath.getSpectralGap();
	sCurrent.dLambda = -1.0;
	sCurrent.nViolations = -1.0;

	printLogLine(LogOutput, sInput, cPath, sCurrent.vdCoefficients, sCurrent.nInfinities);

	return sCurrent;
}

void optimiseAndLog(Path &cPath, std::ofstream &LogOutput, Solution &sCurrent, InputData &sInput) 
{
	cPath.optimise();
	cPath.saveDistances();

	std::cout << "Barriers: " << cPath.getBarriers() << ", Spectral gap: " << cPath.getSpectralGap() << std::endl;

	sCurrent.nBarriers = cPath.getBarriers();
	sCurrent.dSpectralGap = cPath.getSpectralGap();
	sCurrent.dLambda = cPath.getLambda();
	sCurrent.nViolations = cPath.getViolations();
	sCurrent.vnSnapshots = cPath.getSnapshots();

	printLogLine(LogOutput, sInput, cPath, sCurrent.vdCoefficients, sCurrent.nInfinities);
}
