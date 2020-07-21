/*
 * The C++ implementation of the COMet script
 * Also allows to perform SGOOP (the original Tiwary approach, without path)
 * Author: Ladislav Hovan
 */

#include "main.h"

int main(int argc, char *argv[]) 
{
	string strInputFile = loadCommandLineInput(argc, argv);

	InputData sInput = loadInput(strInputFile);
	checkInput(sInput);

	std::cout.precision(sInput.nPrecision - 1);
	std::cout << std::scientific;

	std::cout << "Loading the new collective variables data for the proposed path... ";
	vector2d vvdNewColvarDataPath = load2dDataFromFile(sInput.strNewColvarFilePath);
	vector2d vvdNewColvarsPath = selectColumns(vvdNewColvarDataPath, sInput.vnNewColumnsPath, sInput, true);
	std::cout << "Done" << std::endl;

	// Some vectors are created here but not filled unless needed
	vector2d vvdNewColvars;
	vector<double> vdSnapshotEnergies, vdSummedBias, vdEbetacList;

	std::cout << "Loading the original metadynamics data... ";
	vector2d vvdMetadColvarData = load2dDataFromFile(sInput.strMetadColvarFile);
	std::cout << "Done" << std::endl;

	std::cout << "Loading the new collective variables data... ";
	vector2d vvdNewColvarData = load2dDataFromFile(sInput.strNewColvarFile);
	std::cout << "Done" << std::endl;

	// Make sure that the bias and new collective variables are aligned and have the same number of rows
	alignRows(vvdMetadColvarData, vvdNewColvarData, sInput.nTimeColumn);

	vector<int> vnAllBiasColumns = sInput.vnBiasColumns;
	vnAllBiasColumns.insert(vnAllBiasColumns.end(), sInput.vnRBiasColumns.begin(), sInput.vnRBiasColumns.end());
	vdSummedBias = sumColumns(vvdMetadColvarData, vnAllBiasColumns);
	vvdNewColvars = selectColumns(vvdNewColvarData, sInput.vnNewColumns, sInput, false);

	// Calculate or load work estimates if rbias not included
	if (sInput.vnRBiasColumns.empty())
		vdEbetacList = loadEbetacList(sInput);

	vector2d vvdLimits;
	if (sInput.bRescale)
	{
		std::cout << "Rescaling the collective variables... ";

		// The order is important, otherwise path doesn't get rescaled
		vvdLimits = getLimits(vvdNewColvars);
		vvdNewColvarsPath = rescaleDataRange(vvdNewColvars, vvdNewColvarsPath, sInput, true);
		vvdNewColvars = rescaleDataRange(vvdNewColvars, vvdNewColvars, sInput, false);

		std::cout << "Done" << std::endl;
	}
	
	backupFile(sInput.strLogFile);
	std::ofstream LogOutput(sInput.strLogFile);
	// Not used if not needed
	std::ofstream PeriodicLog;

	if (sInput.bConvertPeriodic)
	{
		backupFile("log_periodic_coeffs.dat");
		PeriodicLog.open("log_periodic_coeffs.dat");
	}
	
	vector<double> vdCoefficients;
	vector<double> vdPeriodicCoeffs;
	vector2d vvdCoefficientsList;  // May not be used
	vector2d vvdPerCoeffList;  // Same

	if (sInput.bTryAll)
	{
		vvdCoefficientsList = createTestingList(sInput.vnNewColumns.size());
		sInput.nAnnealingSteps = vvdCoefficientsList.size();
		vdCoefficients = vvdCoefficientsList[0];
		setSignificant(vdCoefficients, sInput.nPrecision);
	}
	else
	{
		if (sInput.strCoeffList.empty())
		{
			vdCoefficients = initialiseCoefficients(sInput.vnNewColumns.size());
			setSignificant(vdCoefficients, sInput.nPrecision);

			if (sInput.bConvertPeriodic)
				vdPeriodicCoeffs = vector<double>(sInput.vnPeriodicColumns.size(), 0.0);
		}
		else
		{
			vvdCoefficientsList = load2dDataFromFile(sInput.strCoeffList);
			sInput.nAnnealingSteps = vvdCoefficientsList.size();
			vdCoefficients = vvdCoefficientsList[0];
			setSignificant(vdCoefficients, sInput.nPrecision);

			if (sInput.bConvertPeriodic)
			{
				vvdPerCoeffList = load2dDataFromFile(sInput.strPeriodicCoeffList);

				if (vvdPerCoeffList.size() != vvdCoefficientsList.size())
				{
					std::cerr << "The numbers of coefficients and periodic coefficients provided don't match"
						<< std::endl;
					exit (INVALID_INPUT);
				}

				vdPeriodicCoeffs = vvdPerCoeffList[0];
				setSignificant(vdPeriodicCoeffs, sInput.nPrecision);
			}
		}
	}

	Solution sBest;
	sBest.vdCoefficients = vdCoefficients;
	sBest.dSpectralGap = -1.0;

	Solution sPrevious = sBest;

	// Pseudorandom number generator
	std::mt19937 Mersenne(0);

	if (sInput.nSeed < 0)
	{
		std::random_device Random;  // Random seed
		Mersenne.seed(Random());
	}
	else
		Mersenne.seed(sInput.nSeed);

	std::uniform_real_distribution<double> Distribution(0.0, 1.0);

	for (int nStep = 1; nStep <= sInput.nAnnealingSteps; ++nStep, decreasekT(sInput))
	{
		std::cout << "\nAnnealing step " << nStep << ", kT = " << sInput.dAnnealingkT << std::endl;
		std::cout << "The coefficients are: ";
		printVector(vdCoefficients, true);

		Solution sCurrent = sInput.bNoPath ?
			tryCoefficients(vdCoefficients, vdPeriodicCoeffs, vdSummedBias, vvdNewColvars,
				vdEbetacList, sInput, LogOutput) :
			tryCoefficientsPath(vdCoefficients, vdPeriodicCoeffs, vdSummedBias, vvdNewColvars, 
				vvdNewColvarsPath, vdEbetacList, sInput, LogOutput);

		if ((sCurrent.nViolations < sPrevious.nViolations) || ((sCurrent.nViolations == sPrevious.nViolations) &&
			(Distribution(Mersenne) < exp((sCurrent.dSpectralGap - sPrevious.dSpectralGap) / sInput.dAnnealingkT))))
		{
			sPrevious = sCurrent;

			if (sCurrent.dSpectralGap > sBest.dSpectralGap)
				sBest = sCurrent;
		}

		// Log periodic coefficients if used
		if (sInput.bConvertPeriodic)
		{
			for (auto &dCoeff: vdPeriodicCoeffs)
				PeriodicLog << dCoeff << " ";

			PeriodicLog << "\n";
		}

		// Perturb the coefficients from sPrevious
		if (vvdCoefficientsList.empty())
		{
			if (sInput.bRandom)
			{
				vdCoefficients = randomiseCoefficients(vdCoefficients.size(), Mersenne, sInput.vdLowerLimits, 
					sInput.vdUpperLimits);
				setSignificant(vdCoefficients, sInput.nPrecision);

				if (sInput.bConvertPeriodic)
				{
					vdPeriodicCoeffs = randomisePeriodicCoefficients(vdPeriodicCoeffs.size(), Mersenne,
						sInput.vdPeriodicRanges);
					setSignificant(vdPeriodicCoeffs, sInput.nPrecision);
				}
			}
			else
			{
				vdCoefficients = perturbCoefficients(sPrevious.vdCoefficients, Mersenne, sInput.dPerturbationLimit, 
					sInput.vdLowerLimits, sInput.vdUpperLimits);
				setSignificant(vdCoefficients, sInput.nPrecision);

				if (sInput.bConvertPeriodic)
				{
					vdPeriodicCoeffs = perturbPeriodicCoefficients(sPrevious.vdPeriodicCoeffs, Mersenne,
						sInput.vdPeriodicRanges, sInput.dPerturbationLimit);
					setSignificant(vdPeriodicCoeffs, sInput.nPrecision);
				}
			}
		}
		// Or load the next one from the list
		else if (nStep < sInput.nAnnealingSteps)
		{
			vdCoefficients = vvdCoefficientsList[nStep];
			setSignificant(vdCoefficients, sInput.nPrecision);

			if (sInput.bConvertPeriodic)
			{
				vdPeriodicCoeffs = vvdPerCoeffList[nStep];
				setSignificant(vdPeriodicCoeffs, sInput.nPrecision);
			}
		}
	}

	LogOutput.close();
	PeriodicLog.close();

	std::cout << "\nThe best solution:" << std::endl;
	std::cout << "Barriers: " << sBest.nBarriers << ", Spectral gap: " << sBest.dSpectralGap << std::endl;
	std::cout << "Coefficients: ";
	printVector(sBest.vdCoefficients, true);

	if (sInput.bConvertPeriodic)
	{
		std::cout << "Periodic conversion coefficients: ";
		printVector(sBest.vdPeriodicCoeffs, true);
	}

	std::cout << "Lambda: " << sBest.dLambda << std::endl;
	std::cout << "Best snapshots: ";
	printVector(sBest.vnSnapshots, true);

	if (sInput.bStrictOrder) 
		std::cout << "Violations: " << sBest.nViolations << std::endl;

	std::cout << "Infinities in reweighted FES: " << sBest.nInfinities << std::endl;

	if (!sInput.bNoPath)
	{
		// Create a file detailing path progression
		vector<double> vdFirst, vdSecond;
		vdFirst.reserve(vvdNewColvarsPath.size());

		if (sInput.dZLimit >= 0.0 || sInput.bPath2D)
			vdSecond.reserve(vvdNewColvarsPath.size());

		pvd pvdPathValues(vdFirst, vdSecond);
		vector2d vvdConvertedPath;

		if (sInput.bConvertPeriodic)
		{
			vvdConvertedPath = convertPeriodic(vvdNewColvarsPath, sBest.vdPeriodicCoeffs, sInput);

			calculatePath(vvdConvertedPath, vvdConvertedPath, sBest.vdCoefficients, sBest.vnSnapshots, sBest.dLambda, 
				sInput, pvdPathValues);
		}
		else
			calculatePath(vvdNewColvarsPath, vvdNewColvarsPath, sBest.vdCoefficients, sBest.vnSnapshots, sBest.dLambda,
				sInput, pvdPathValues);

		backupFile("s_progression.dat");
		writeDataFile("s_progression.dat", pvdPathValues.first);

		// Create a file with PCV values for the original simulation
		pvdPathValues.first.clear();
		pvdPathValues.first.reserve(vvdNewColvars.size());
		pvdPathValues.second.clear();
		pvdPathValues.second.reserve(vvdNewColvars.size());

		sInput.bPath2D = true;
		vector2d vvdConverted;

		if (sInput.bConvertPeriodic)
		{
			vvdConverted = convertPeriodic(vvdNewColvars, sBest.vdPeriodicCoeffs, sInput);

			calculatePath(vvdConvertedPath, vvdConverted, sBest.vdCoefficients, sBest.vnSnapshots, sBest.dLambda,
				sInput, pvdPathValues);
		}
		else
			calculatePath(vvdNewColvarsPath, vvdNewColvars, sBest.vdCoefficients, sBest.vnSnapshots, sBest.dLambda,
				sInput, pvdPathValues);

		vector2d vvdPathValues;
		for (int nPos = 0; nPos < vvdNewColvars.size(); ++nPos)
			vvdPathValues.push_back({ pvdPathValues.first[nPos], pvdPathValues.second[nPos] });

		backupFile("pcv_values.dat");
		writeDataFile("pcv_values.dat", vvdPathValues);

		// Create a Colvar file
		if (sInput.bDescale)
		{
			vvdNewColvarDataPath = load2dDataFromFile(sInput.strNewColvarFilePath);
			vvdNewColvarsPath = selectColumns(vvdNewColvarDataPath, sInput.vnNewColumnsPath, sInput, true);

			if (sInput.bConvertPeriodic)
			{
				vdPeriodicCoeffs = sBest.vdPeriodicCoeffs;

				for (int nCV = 0; nCV < sBest.vdCoefficients.size(); ++nCV)
					vdPeriodicCoeffs[nCV] *= (vvdLimits[1][nCV] - vvdLimits[0][nCV]);

				vector2d vvdConvertedPath = convertPeriodic(vvdNewColvarsPath, vdPeriodicCoeffs, sInput);
				createColvarFile(vvdConvertedPath, sBest.vnSnapshots, sInput.nPrecision);
			}
			else
				createColvarFile(vvdNewColvarsPath, sBest.vnSnapshots, sInput.nPrecision);
		}
		else
		{
			if (sInput.bConvertPeriodic)
			{
				vector2d vvdConvertedPath = convertPeriodic(vvdNewColvarsPath, sBest.vdPeriodicCoeffs, sInput);
				createColvarFile(vvdConvertedPath, sBest.vnSnapshots, sInput.nPrecision);
			}
			else
				createColvarFile(vvdNewColvarsPath, sBest.vnSnapshots, sInput.nPrecision);
		}
	}
	
	// Create a Plumed file
	createPlumedFile(sBest, vvdLimits, sInput);

	return NO_ERROR;
}
