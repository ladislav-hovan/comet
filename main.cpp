/*
 * The C++ implementation of the COMet script
 * Also allows to perform SGOOP (the original Tiwary approach, without path)
 * Author: Ladislav Hovan
 */

#include "main.h"

int main(int argc, char *argv[]) 
{
	std::string strInputFile = loadCommandLineInput(argc, argv);

	InputData sInput = loadInput(strInputFile);
	checkInput(sInput);

	std::cout << "Loading the new collective variables data for the proposed path... ";
	vector2d vvdNewColvarDataPath = load2dDataFromFile(sInput.strNewColvarFilePath);
	vector2d vvdNewColvarsPath = selectColumns(vvdNewColvarDataPath, sInput.vnNewColumnsPath, sInput, true);
	std::cout << "Done" << std::endl;

	// Some vectors are created here but not filled unless needed
	vector2d vvdNewColvars;
	vector<double> vdSnapshotEnergies;
	vector<double> vdSummedBias;
	vector<double> vdEbetacList;
	vector<double> vdScalings;

	std::cout << "Loading the original metadynamics data... ";
	vector2d vvdMetadColvarData = load2dDataFromFile(sInput.strMetadColvarFile);
	std::cout << "Done" << std::endl;

	std::cout << "Loading the new collective variables data... ";
	vector2d vvdNewColvarData = load2dDataFromFile(sInput.strNewColvarFile);
	std::cout << "Done" << std::endl;

	// Make sure that the bias and new collective variables are aligned and have the same number of rows
	alignRows(vvdMetadColvarData, vvdNewColvarData, sInput.nTimeColumn);

	vdSummedBias = sumColumns(vvdMetadColvarData, sInput.vnBiasColumns);
	vvdNewColvars = selectColumns(vvdNewColvarData, sInput.vnNewColumns, sInput, false);

	// Either get work estimates or take weights of the hills which incorporate them
	if (sInput.nScalingColumn == -1)
		vdEbetacList = loadEbetacList(sInput);
	else
	{
		for (auto &row: vvdMetadColvarData)
			vdScalings.push_back(row[sInput.nScalingColumn]);
	}

	if (sInput.bRescale)
	{
		std::cout << "Rescaling the collective variables... ";

		// The order is important, otherwise path doesn't get rescaled
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
	}
	else
	{
		if (sInput.strCoeffList.empty())
		{
			vdCoefficients = initialiseCoefficients(sInput.vnNewColumns.size());

			if (sInput.bConvertPeriodic)
				vdPeriodicCoeffs = vector<double>(sInput.vnPeriodicColumns.size(), 0.0);
		}
		else
		{
			vvdCoefficientsList = load2dDataFromFile(sInput.strCoeffList);
			sInput.nAnnealingSteps = vvdCoefficientsList.size();
			vdCoefficients = vvdCoefficientsList[0];

			if (sInput.bConvertPeriodic)
			{
				vvdPerCoeffList = load2dDataFromFile(sInput.strPeriodicCoeffList);

				if (vvdPerCoeffList.size() != vvdCoefficientsList.size())
				{
					std::cerr << "The numbers of coefficients and periodic coefficients provided don't match" << std::endl;
					exit (INVALID_INPUT);
				}

				vdPeriodicCoeffs = vvdPerCoeffList[0];
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
				vdEbetacList, vdScalings, sInput, LogOutput) :
			tryCoefficientsPath(vdCoefficients, vdPeriodicCoeffs, vdSummedBias, vvdNewColvars, 
				vvdNewColvarsPath, vdEbetacList, vdScalings, sInput, LogOutput);

		if ((sCurrent.nViolations < sPrevious.nViolations) || ((sCurrent.nViolations == sPrevious.nViolations) &&
			(Distribution(Mersenne) < exp((sCurrent.dSpectralGap - sPrevious.dSpectralGap) / sInput.dAnnealingkT))))
		{
			sPrevious = sCurrent;

			if (sCurrent.dSpectralGap > sBest.dSpectralGap)
			{
				sBest = sCurrent;

				if (fs::exists("dist_matrix.dat"))
					fs::rename("dist_matrix.dat", "best_dist_matrix.dat");
			}
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
				vdCoefficients = randomiseCoefficients(vdCoefficients.size(), Mersenne);

				if (sInput.bConvertPeriodic)
					vdPeriodicCoeffs = randomisePeriodicCoefficients(vdPeriodicCoeffs.size(), Mersenne, sInput.vdPeriodicRanges);
			}
			else
			{
				vdCoefficients = perturbCoefficients(sPrevious.vdCoefficients, Mersenne);

				if (sInput.bConvertPeriodic)
					vdPeriodicCoeffs = perturbPeriodicCoefficients(sPrevious.vdPeriodicCoeffs, Mersenne, sInput.vdPeriodicRanges);
			}
		}
		// Or load the next one from the list
		else if (nStep < sInput.nAnnealingSteps)
		{
			vdCoefficients = vvdCoefficientsList[nStep];

			if (sInput.bConvertPeriodic)
				vdPeriodicCoeffs = vvdPerCoeffList[nStep];
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

	if (sInput.bConvertPeriodic)
	{
		vector2d vvdConvertedPath = convertPeriodic(vvdNewColvarsPath, sBest.vdPeriodicCoeffs, sInput);
		createColvarFile(vvdConvertedPath, sBest.vnSnapshots);
	}
	else
		createColvarFile(vvdNewColvarsPath, sBest.vnSnapshots);

	return NO_ERROR;
}
