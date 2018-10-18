/*
 * path.cpp
 *
 *  Created on: 17 Mar 2016
 *      Author: ladislav
 */

#include "path.h"

Path::Path() 
{
	// Auto-generated constructor stub
}

Path::Path(vector<double> &vdFes, InputData &sInput, double dBestEnergy) 
{
	// Special constructor for the s variable, no optimisation needed, just SG calculation

	m_vdEnergies = vdFes;
	m_nTotalSnapshots = vdFes.size();
	m_nLengthTarget = vdFes.size();
	m_nLengthTolerance = 0;
	m_nBarriersSet = sInput.nBarriers;
	m_dkT = sInput.dkT;
	m_dThreshold = sInput.dThreshold * sInput.dkT;
	m_bSmoothCount = sInput.bSmoothCount;
	sBestPath.dEnergy = dBestEnergy;
	initialiseSnapshots(m_nLengthTarget, false);
}

Path::Path(vector<double> &vdPathEnergies, vector<double> &vdCoefficients, InputData &sInput,
	vector2d &vvdNewColvarsPath): Path::Path(vdCoefficients, sInput, vvdNewColvarsPath) 
{
	// Constructor to be used for direct assignment of energies

	m_vdEnergies = vdPathEnergies;
}

Path::Path(vector<double> &vdGrid, vector<double> &vdFes, vector<double> &vdCoefficients, InputData &sInput,
		   vector2d &vvdNewColvarsPath, bool bAssign): Path::Path(vdCoefficients, sInput, vvdNewColvarsPath) 
{
	// Constructor to be used for energies from reweighted FES

	if (bAssign)
		assignEnergies(vdGrid, vdFes, vvdNewColvarsPath);
}

Path::Path(vector<double> &vdCoefficients, InputData &sInput, vector2d &vvdNewColvarsPath) 
{
	// The private constructor to be used by other constructors

	m_nLengthTarget = sInput.nLengthTarget;
	m_nLengthTolerance = sInput.nLengthTolerance;
	m_nTotalSnapshots = vvdNewColvarsPath.size();
	m_nBarriersSet = sInput.nBarriers;
	m_dkT = sInput.dkT;
	m_dThreshold = sInput.dThreshold * sInput.dkT;
	m_bStrict = sInput.bStrictOrder;
	m_bSmoothCount = sInput.bSmoothCount;
	m_vdCoefficients = vdCoefficients;
	setPeriodicity(sInput);
	loadMatrix(vdCoefficients, vvdNewColvarsPath, sInput.bEuclidean);
	initialiseSnapshots(m_nLengthTarget, false);  // optimise() will look at other lengths too
}

Path::~Path() 
{
	// Auto-generated destructor stub
}

void Path::initialiseSnapshots(unsigned int nSnapshots, bool bPrint) 
{
	m_vnSnapshots.clear();
	m_vnSnapshots.push_back(0);

	for (unsigned int nCount = 1; nCount < nSnapshots - 1; ++nCount)
		// Integer division is intentional here
		m_vnSnapshots.push_back((nCount * m_nTotalSnapshots) / (nSnapshots - 1));

	m_vnSnapshots.push_back(m_nTotalSnapshots - 1);

	if (bPrint)
	{
		std::cout << "Initial guess at a path: ";
		printPath();
	}
}

void Path::printPath() 
{
	for (unsigned int nCount = 0; nCount < m_vnSnapshots.size(); ++nCount)
	{
		std::cout << m_vnSnapshots.at(nCount);

		if (nCount < m_vnSnapshots.size() - 1)
			std::cout << ", ";
	}

	std::cout << std::endl;
}

void Path::setPeriodicity(InputData &sInput) 
{
	m_vdPeriodicRanges.clear();

	for (unsigned int nColumn = 0; nColumn < sInput.vnNewColumnsPath.size(); ++nColumn)
	{
		if (sInput.bConvertPeriodic)
		{
			// Periodic variables were converted, so none remain
			for (unsigned int nCount = 0; nCount < sInput.vnPeriodicColumnsPath.size(); ++nCount)
				m_vdPeriodicRanges.push_back(-1.0);
		}
		else
		{
			int nPeriodic = -1;

			for (unsigned int nCount = 0; nCount < sInput.vnPeriodicColumnsPath.size(); ++nCount)
				// Cast to prevent signed/unsigned comparison
				if (static_cast<signed int>(nColumn) == sInput.vnPeriodicColumnsPath[nCount])
					nPeriodic = nCount;

			if (nPeriodic < 0)
				m_vdPeriodicRanges.push_back(-1.0);
			else
				m_vdPeriodicRanges.push_back(sInput.vdPeriodicRanges[nPeriodic]);
		}
	}
}

void Path::loadMatrix(vector<double> &vdCoefficients, vector2d &vvdNewColvarsPath, bool bEuclidean) 
{
	std::cout << "Loading a matrix of distance values... ";

	int nSize = vvdNewColvarsPath.size();
	vector<double> vdTempRow(nSize, 0.0);

	for (int nCount = 0; nCount < nSize; ++nCount)
		m_vvdMatrix.push_back(vdTempRow);

	// Speed is important, so I won't use bound checking via at()
	for (int nRow = 0; nRow < nSize; ++nRow)
	{
		// Only the two nearest neighbours calculated because only those distances are used
		for (int nColumn = nRow; nColumn < nSize; ++nColumn)
		{
			if (nRow == nColumn)
				m_vvdMatrix[nRow][nColumn] = 0.0;
			else
			{
				double dTempValue = 0.0;

				if (bEuclidean)
				{
					for (unsigned int nColvar = 0; nColvar < vdCoefficients.size(); ++nColvar)
					{
						double dDiff = vvdNewColvarsPath[nRow][nColvar] - vvdNewColvarsPath[nColumn][nColvar];

						if (m_vdPeriodicRanges[nColvar] > 0)
						{
							if (dDiff > m_vdPeriodicRanges[nColvar] / 2)
								dDiff -= m_vdPeriodicRanges[nColvar];
							else if (dDiff < - m_vdPeriodicRanges[nColvar] / 2)
								dDiff += m_vdPeriodicRanges[nColvar];
						}

						dTempValue += pow(vdCoefficients[nColvar] * dDiff, 2);
					}

					dTempValue = sqrt(dTempValue);
				}
				else
				{
					for (unsigned int nColvar = 0; nColvar < vdCoefficients.size(); ++nColvar)
					{
						double dDiff = vvdNewColvarsPath[nRow][nColvar] - vvdNewColvarsPath[nColumn][nColvar];

						if (m_vdPeriodicRanges[nColvar] > 0)
						{
							if (dDiff > m_vdPeriodicRanges[nColvar] / 2)
								dDiff -= m_vdPeriodicRanges[nColvar];
							else if (dDiff < - m_vdPeriodicRanges[nColvar] / 2)
								dDiff += m_vdPeriodicRanges[nColvar];
						}

						dTempValue += fabs(vdCoefficients[nColvar] * dDiff);
					}
				}

				m_vvdMatrix[nRow][nColumn] = dTempValue;
				m_vvdMatrix[nColumn][nRow] = dTempValue;
			}
		}
	}

	std::cout << "Done" << std::endl;
}

void Path::assignEnergies(vector<double> &vdGrid, vector<double> &vdFes, vector2d &vvdNewColvarsPath) 
{
	m_vdEnergies.reserve(vvdNewColvarsPath.size());

	for (unsigned int nSnapshot = 0; nSnapshot < vvdNewColvarsPath.size(); ++nSnapshot)
	{
		double dTotalValue = 0.0;

		for (unsigned int nColumn = 0; nColumn < vvdNewColvarsPath[nSnapshot].size(); ++nColumn)
			dTotalValue += vvdNewColvarsPath[nSnapshot][nColumn] * m_vdCoefficients[nColumn];

		if (dTotalValue >= vdGrid[0] && dTotalValue <= vdGrid[vdGrid.size() - 1])
		{
			for (unsigned int nValue = 0; nValue < vdGrid.size() - 1; ++nValue)
			{
				if (vdGrid[nValue] <= dTotalValue && vdGrid[nValue + 1] >= dTotalValue)
				{
					double dNumerator = vdFes[nValue] * (dTotalValue - vdGrid[nValue]) + 
										vdFes[nValue + 1] * (vdGrid[nValue + 1] - dTotalValue);
					double dDenominator = vdGrid[nValue + 1] - vdGrid[nValue];
					m_vdEnergies.push_back(dNumerator / dDenominator);
					break;
				}
			}
		}
		else
		{
			std::cerr << "The value of collective variable for path is outside the grid range" << std::endl;
			exit (OUTSIDE_GRID);
		}
	}

	writeDataFile("reweighted_fes.dat", m_vdEnergies);
}

void Path::assignEnergies(vector2d &vvdGridList, vector2d &vvdFesList, vector2d &vvdNewColvarsPath) 
{
	m_vdEnergies.reserve(vvdNewColvarsPath.size());

	for (unsigned int nSnapshot = 0; nSnapshot < vvdNewColvarsPath.size(); ++nSnapshot)
	{
		vector<double> vdEnergyValues;

		for (unsigned int nColumn = 0; nColumn < vvdNewColvarsPath[nSnapshot].size(); ++nColumn)
		{
			double dColvarValue = vvdNewColvarsPath[nSnapshot][nColumn];

			if (dColvarValue >= vvdGridList[nColumn][0] && dColvarValue <= vvdGridList[nColumn][vvdGridList[nColumn].size() - 1])
			{
				for (unsigned int nValue = 0; nValue < vvdGridList[nColumn].size() - 1; ++nValue)
				{
					if (vvdGridList[nColumn][nValue] <= dColvarValue && vvdGridList[nColumn][nValue + 1] >= dColvarValue)
					{
						double dNumerator = vvdFesList[nColumn][nValue] * (dColvarValue - vvdGridList[nColumn][nValue]) +
											vvdFesList[nColumn][nValue + 1] * (vvdGridList[nColumn][nValue + 1] - dColvarValue);
						double dDenominator = vvdGridList[nColumn][nValue + 1] - vvdGridList[nColumn][nValue];
						vdEnergyValues.push_back(dNumerator / dDenominator);
						break;
					}
				}
			}
			else
			{
				std::cerr << "The value of collective variable for path is outside the grid range" << std::endl;
				exit (OUTSIDE_GRID);
			}
		}

		double dTotalValue = 0.0;

		for (unsigned int nColvar = 0; nColvar < m_vdCoefficients.size(); ++nColvar)
			dTotalValue += m_vdCoefficients[nColvar] * vdEnergyValues[nColvar];

		m_vdEnergies.push_back(dTotalValue);
	}

	writeDataFile("reweighted_fes.dat", m_vdEnergies);
}

void Path::optimise() 
{
	for (unsigned int nLength = m_nLengthTarget - m_nLengthTolerance; nLength <= m_nLengthTarget + m_nLengthTolerance; ++nLength)
	{
		std::cout << ">>>\nTesting path of length " << nLength << std::endl;
		initialiseSnapshots(nLength);
		optimiseSetLength();

		double dEnergy = getEnergy();

		// If strict order is not enforced, all violations are zero so the middle bit is never true
		if (sBestPath.vnSnapshots.empty() || m_nViolations < sBestPath.nViolations ||
			(dEnergy < sBestPath.dEnergy && m_nViolations == sBestPath.nViolations))
		{
			sBestPath.vnSnapshots = m_vnSnapshots;
			sBestPath.dLambda = 2.3 / pow(m_dFirstAverage, 2);
			sBestPath.dEnergy = dEnergy;
			sBestPath.nViolations = m_nViolations;
		}
	}

	m_vnSnapshots = sBestPath.vnSnapshots;
	std::cout << ">>>\nUsing " << m_vnSnapshots.size() << " snapshots" << std::endl;
}

int Path::checkAllDistances() 
{
	int nViolations = 0;

	for (unsigned int nPosition = 0; nPosition < m_vnSnapshots.size(); ++nPosition)
		nViolations += checkDistances(nPosition);

	return nViolations;
}

int Path::checkDistances(unsigned int nPosition) 
{
	int nViolations = 0;

	// Get the first second distance and compare with relevant first ones
	if (nPosition > 1)
	{
		double dSecondDist = m_vdSecondDistances[nPosition - 2];

		if (m_vdFirstDistances[nPosition - 2] > dSecondDist)
			++nViolations;
		if (m_vdFirstDistances[nPosition - 1] > dSecondDist) 
			++nViolations;
	}

	// Get the other one and compare as before
	if (nPosition < m_vnSnapshots.size() - 2)
	{
		double dSecondDist = m_vdSecondDistances[nPosition];

		if (m_vdFirstDistances[nPosition] > dSecondDist) 
			++nViolations;
		if (m_vdFirstDistances[nPosition + 1] > dSecondDist) 
			++nViolations;
	}

	return nViolations;
}

void Path::optimiseSetLength() 
{
	initialiseOptimisationVars();

	// At max 100 steps of optimisation
	for (int nStep = 1; nStep <= 100; ++nStep)
	{
		vector<int> vnBackupSnapshots = m_vnSnapshots;
		int nTotalViolations = m_bStrict ? checkAllDistances() : 0;

		for (unsigned int nPosition = 1; nPosition < m_vnSnapshots.size() - 1; ++nPosition)
		{
			double dBestEnergy = dInf;
			unsigned int nBestSnapshot = m_vnSnapshots[nPosition - 1] + 1;

			for (int nSnapshot = m_vnSnapshots[nPosition - 1] + 1; nSnapshot < m_vnSnapshots[nPosition + 1]; ++nSnapshot)
			{
				changeAndRecalculate(nPosition, nSnapshot);
				double dEnergy = getEnergy();
				bool bAccept = false;

				if (m_bStrict)
				{
					int nCurrentViolations = checkAllDistances();

					// This way, the number of violations can never increase
					if (nCurrentViolations > nTotalViolations)
						continue;
					// This ensures it will decrease if it can
					else if (nCurrentViolations < nTotalViolations)
					{
						bAccept = true;
						nTotalViolations = nCurrentViolations;
					}
				}

				if (bAccept || dEnergy < dBestEnergy)
				{
					dBestEnergy = dEnergy;
					nBestSnapshot = nSnapshot;
				}
			}

			changeAndRecalculate(nPosition, nBestSnapshot);
		}

		if (vnBackupSnapshots == m_vnSnapshots)
		{
			m_nViolations = nTotalViolations;
			std::cout << "Deterministic minimisation yields no further improvement after " << nStep << " steps" << std::endl;
			std::cout << "The final energy is " << getEnergy() << std::endl;
			return;
		}
		else
		{
			m_nViolations = nTotalViolations;
			std::cout << "New best path: ";
			printPath();

			if (m_bStrict) 
				std::cout << "Violations: " << nTotalViolations << std::endl;
		}
	}

	std::cout << "Deterministic minimisation has not converged after 100 steps, taking current result" << std::endl;
	std::cout << "The final energy is " << getEnergy() << std::endl;
}

double Path::getEnergy() 
{
	double dFirstDev = calculateDeviation(m_vdFirstDistances, m_dFirstAverage);
	double dSecondDev = calculateDeviation(m_vdSecondDistances, m_dSecondAverage);

	// An empirical energy function for optimisation
	return log(dFirstDev + 0.35 * m_dFirstAverage + 0.65 * dSecondDev + 0.1 * m_dSecondAverage);
}

double Path::getBestEnergy() 
{
	return sBestPath.dEnergy;
}

double Path::recomputeEnergy() 
{
	initialiseOptimisationVars();

	return getEnergy();
}

double Path::getSpectralGap() 
{
	if (m_dSpectralGap < 0.0)
		calculateSpectralGap();

	return m_dSpectralGap;
}

int Path::getBarriers()
{
	if (m_dSpectralGap < 0.0)
		calculateSpectralGap();

	return m_nBarriers;
}

int Path::getViolations()
{
	if (sBestPath.vnSnapshots.empty())
		return 0;
	else
		return sBestPath.nViolations;
}

double Path::getLambda() 
{
	if (sBestPath.vnSnapshots.empty())
	{
		initialiseOptimisationVars();

		return 2.3 / pow(m_dFirstAverage, 2);
	}
	else
		return sBestPath.dLambda;
}

vector<double> Path::getEigenvalues(int nEigenvalues) 
{
	if (m_dSpectralGap < 0.0)
		calculateSpectralGap();

	// The cast is only there to stop warnings
	if (nEigenvalues < 0 || nEigenvalues >= static_cast<signed>(m_vdEigenvalues.size()))
		return m_vdEigenvalues;
	else
	{
		vector<double> vdSubset(m_vdEigenvalues.begin(), m_vdEigenvalues.begin() + nEigenvalues);

		return vdSubset;
	}
}

vector<int> Path::getSnapshots() 
{
	if (sBestPath.vnSnapshots.empty())
	{
		initialiseSnapshots(m_nLengthTarget, false);
		return m_vnSnapshots;
	}
	else
		return sBestPath.vnSnapshots;
}

void Path::saveDistances(std::string strFilename) 
{
	vector2d vvdMatrix;

	for (int &nRow: m_vnSnapshots)
	{
		vector<double> vdRow;

		for (int &nCol: m_vnSnapshots)
			vdRow.push_back(m_vvdMatrix[nRow][nCol]);

		vvdMatrix.push_back(vdRow);
	}

	writeDataFile(strFilename, vvdMatrix);
}

void Path::initialiseOptimisationVars() 
{
	m_vdFirstDistances.clear();
	m_vdSecondDistances.clear();

	double dFirstSum = 0.0;
	double dSecondSum = 0.0;

	for (unsigned int nCount = 0; nCount < m_vnSnapshots.size() - 1; ++nCount)
	{
		double dFirstDist = m_vvdMatrix[m_vnSnapshots[nCount]][m_vnSnapshots[nCount + 1]];
		dFirstSum += dFirstDist;
		m_vdFirstDistances.push_back(dFirstDist);

		if (nCount < m_vnSnapshots.size() - 2)
		{
			double dSecondDist = m_vvdMatrix[m_vnSnapshots[nCount]][m_vnSnapshots[nCount + 2]];
			dSecondSum += dSecondDist;
			m_vdSecondDistances.push_back(dSecondDist);
		}
	}

	m_dFirstAverage = dFirstSum / m_vdFirstDistances.size();
	m_dSecondAverage = dSecondSum / m_vdSecondDistances.size();
}

void Path::changeAndRecalculate(unsigned int nPosition, int nSnapshot) 
{
	if (nPosition > 0)
	{
		double dNewFirst = m_vvdMatrix[m_vnSnapshots[nPosition - 1]][nSnapshot];
		m_dFirstAverage += (dNewFirst - m_vdFirstDistances[nPosition - 1]) / m_vdFirstDistances.size();
		m_vdFirstDistances[nPosition - 1] = dNewFirst;
	}

	if (nPosition > 1)
	{
		double dNewSecond = m_vvdMatrix[m_vnSnapshots[nPosition - 2]][nSnapshot];
		m_dSecondAverage += (dNewSecond - m_vdSecondDistances[nPosition - 2]) / m_vdSecondDistances.size();
		m_vdSecondDistances[nPosition - 2] = dNewSecond;
	}

	if (nPosition < m_vnSnapshots.size() - 1)
	{
		double dNewFirst = m_vvdMatrix[nSnapshot][m_vnSnapshots[nPosition + 1]];
		m_dFirstAverage += (dNewFirst - m_vdFirstDistances[nPosition]) / m_vdFirstDistances.size();
		m_vdFirstDistances[nPosition] = dNewFirst;
	}

	if (nPosition < m_vnSnapshots.size() - 2)
	{
		double dNewSecond = m_vvdMatrix[nSnapshot][m_vnSnapshots[nPosition + 2]];
		m_dSecondAverage += (dNewSecond - m_vdSecondDistances[nPosition]) / m_vdSecondDistances.size();
		m_vdSecondDistances[nPosition] = dNewSecond;
	}

	m_vnSnapshots[nPosition] = nSnapshot;
	m_dSpectralGap = -1.0;
	m_nBarriers = -1;
}

double Path::calculateDeviation(vector<double> &vdData, double dAverage) 
{
	// This assumes average is correct and doesn't recalculate it
	double dVariance = 0.0;

	for (auto dValue: vdData)
		dVariance += pow((dValue - dAverage), 2);

	dVariance /= vdData.size();

	return sqrt(dVariance);
}

double Path::calculateDeviation(vector<double> &vdData) 
{
	double dAverage = 0.0;

	for (auto dValue: vdData)
		dAverage += dValue;

	dAverage /= vdData.size();

	return calculateDeviation(vdData, dAverage);
}

void Path::determineSpectralGap() 
{
	if (m_nBarriersSet >= 0)
	{
		// Using the set number of barriers
		m_nBarriers = m_nBarriersSet;
		m_dSpectralGap = m_vdEigenvalues[m_nBarriers] - m_vdEigenvalues[m_nBarriers + 1];
	}
	else if (m_nBarriersSet == -1)
	{
		// Using the "first real" spectral gap definition
		double dPreviousGap = 0.0;

		for (unsigned int nPos = 0; nPos < m_vdEigenvalues.size() - 2; ++nPos)
		{
			double dCurrentGap = m_vdEigenvalues[nPos] - m_vdEigenvalues[nPos + 1];
			double dNextGap = m_vdEigenvalues[nPos + 1] - m_vdEigenvalues[nPos + 2];

			if (dCurrentGap > dPreviousGap && dCurrentGap > dNextGap)
			{
				m_dSpectralGap = dCurrentGap;
				m_nBarriers = nPos;
				break;
			}
			else
				dPreviousGap = dCurrentGap;
		}
	}
	else if (m_nBarriersSet == -2)
	{
		// Actually trying to count barriers higher than threshold
		if (m_bSmoothCount)
		{
			vector<double> vdSmoothed = smoothSurface(m_vdEnergies, true);
			m_nBarriers = countBarriers(vdSmoothed);
		}
		else
			m_nBarriers = countBarriers(m_vdEnergies);

		m_dSpectralGap = m_vdEigenvalues[m_nBarriers] - m_vdEigenvalues[m_nBarriers + 1];
	}
	else
	{
		std::cerr << "Invalid setting for number of barriers" << std::endl;
		exit (NOT_IMPLEMENTED);
	}
}

int Path::countBarriers(vector<double> &vdFes) 
{
	// Find all the minima, no matter how shallow
	vector<int> vnMinPositions;

	for (unsigned int nPos = 0; nPos < vdFes.size(); ++nPos)
	{
		bool bMin = true;
		double dEnergy = vdFes[nPos];

		if (nPos > 0)
		{
			if (vdFes[nPos - 1] < dEnergy)
			{
				bMin = false;
				continue;
			}
		}

		if (nPos < (vdFes.size() - 1))
		{
			if (vdFes[nPos + 1] < dEnergy)
				bMin = false;
		}

		if (bMin)
			vnMinPositions.push_back(nPos);
	}

	printVector(vnMinPositions, true);

	// Only one minimum (or possibly zero) - no barriers
	if (vnMinPositions.size() <= 1)
	{
		std::cout << "Found no barriers!" << std::endl;
		return 0;
	}

	// Now determine which of these are separated by less than threshold
	int nBarriers = 0;
	unsigned int nMin1 = 0, nMin2 = 1;

	do
	{
		// Some room for optimisation here, the dMax doesn't need to be determined for the whole range again
		double dMax = -dInf;

		for (int nPos = vnMinPositions[nMin1]; nPos < vnMinPositions[nMin2]; ++nPos)
		{
			if (vdFes[nPos] > dMax)
				dMax = vdFes[nPos];
		}

		if (((dMax - vdFes[vnMinPositions[nMin1]]) > m_dThreshold) && ((dMax - vdFes[vnMinPositions[nMin2]]) > m_dThreshold))
		{
			++nBarriers;
			nMin1 = nMin2;
		}
		else if (vdFes[vnMinPositions[nMin1]] > vdFes[vnMinPositions[nMin2]])
			++nMin1;

		++nMin2;
	}
	while (nMin2 < vnMinPositions.size());

	if (nBarriers == 0)
		std::cout << "Found no barriers!" << std::endl;
	else if (nBarriers == 1)
		std::cout << "Found 1 barrier!" << std::endl;
	else
		std::cout << "Found " << nBarriers << " barriers!" << std::endl;

	return nBarriers;
}

vector<double> Path::smoothSurface(vector<double> &vdFes, bool bGaussian) 
{
	writeDataFile("rough_energies.dat", vdFes);
	vector<double> vdSmoothed;

	if (!bGaussian)
	{
		// Use rolling average
		int nWindow = vdFes.size() / 20;
		int nIncluded = nWindow;
		double dSum = 0.0;

		if (nWindow == 0)
		{
			std::cerr << "The number of points in the free energy surface is too small to use rolling average" << std::endl;
			exit (INVALID_INPUT);
		}

		for (unsigned int nPos = 0; nPos < vdFes.size(); ++nPos)
		{
			if (nPos == 0)
			{
				for (int nCount = 0; nCount < nWindow; ++nCount)
					dSum += vdFes[nCount];
			}
			else
			{
				if (nPos > nWindow)
				{
					dSum -= vdFes[nPos - nWindow - 1];
					--nIncluded;
				}
				if (nPos + nWindow < vdFes.size())
				{
					dSum += vdFes[nPos + nWindow];
					++nIncluded;
				}
			}

			vdSmoothed.push_back(dSum / nIncluded);
		}
	}
	else
	{
		// Use gaussian kernel
		double dStd = vdFes.size() * 0.03;

		for (unsigned int nPos = 0; nPos < vdFes.size(); ++nPos)
		{
			double dNumerator = 0.0, dDenominator = 0.0;

			for (unsigned int nCount = 0; nCount < vdFes.size(); ++nCount)
			{
				double dExp = exp(-pow((nPos - nCount) / dStd, 2) / 2);
				dNumerator += dExp * vdFes[nCount];
				dDenominator += dExp;
			}

			vdSmoothed.push_back(dNumerator / dDenominator);
		}
	}

	writeDataFile("smoothed_energies.dat", vdSmoothed);

	return vdSmoothed;
}

void Path::logPathEnergies() 
{
	std::ofstream PathOutput("path_energies.dat");

	for (const auto &nSnapshot: m_vnSnapshots)
		PathOutput << m_vdEnergies[nSnapshot] << "\n";

	PathOutput.close();
}

void Path::populatekMatrix() 
{
	m_kMatrix.resize(m_vnSnapshots.size(), m_vnSnapshots.size());

	for (unsigned int nRow = 0; nRow < m_vnSnapshots.size(); ++nRow)
	{
		for (unsigned int nColumn = 0; nColumn < m_vnSnapshots.size(); ++nColumn)
		{
			if (nColumn == nRow + 1)
			{
				double dValue = sqrt(exp(-(m_vdEnergies[m_vnSnapshots[nColumn]] - m_vdEnergies[m_vnSnapshots[nRow]]) / m_dkT));

				m_kMatrix(nRow, nColumn) = dValue;
				m_kMatrix(nColumn, nRow) = 1.0 / dValue;
			}
			else if (nColumn != nRow - 1)
				m_kMatrix(nRow, nColumn) = 0.0;
		}

		double dSum = 0.0;

		if (nRow > 0)
			dSum += m_kMatrix(nRow, nRow - 1);
		if (nRow < m_vnSnapshots.size() - 1)
			dSum += m_kMatrix(nRow, nRow + 1);

		m_kMatrix(nRow, nRow) = 1.0 - dSum;
	}
}

void printLogLine(std::ofstream &Output, InputData &sInput, Path &cPath, vector<double> &vdCoefficients, int nInfCount) 
{
	if (sInput.bTryAll)
	{
		static int nColvar = 0;
		printNumberToFile(Output, nColvar++);
	}
	else
		printNumbersToFile(Output, vdCoefficients);

	printNumberToFile(Output, cPath.getBarriers());
	printNumberToFile(Output, cPath.getSpectralGap());
	printNumberToFile(Output, nInfCount);

	if (sInput.bLogEnergy)
		printNumberToFile(Output, cPath.getBestEnergy());

	if (sInput.nLogEigenvalues > 0)
	{
		vector<double> vdEigenvalues = cPath.getEigenvalues(sInput.nLogEigenvalues);
		printNumbersToFile(Output, vdEigenvalues);
	}

	Output << "\n";
}
