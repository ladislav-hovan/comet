/*
 * input.cpp
 *
 *  Created on: 10 Mar 2016
 *      Author: ladislav
 */

#include "input.h"

std::string loadCommandLineInput(int argc, char *argv[]) 
{
	std::string strInputFile;

	if (argc < 2)
		strInputFile = "parameters.input";
	else
	{
		strInputFile = argv[1];
		if (argc > 2)
			std::cout << "More than one command line argument provided, ignoring the rest" << std::endl;
	}

	return strInputFile;
}

vector<double> loadFesFiles(InputData &sInput) 
{
	vector<double> vdEbetacList;

	for (int nFile = 0; nFile < sInput.nFesFilesCount; ++nFile)
	{
		std::stringstream Stream;
		Stream << sInput.strFesPrefix << nFile << ".dat";
		std::string strFilename = Stream.str();
		double dEbetac = getEbetac(strFilename, sInput);
		vdEbetacList.push_back(dEbetac);
	}

	backupFile("generated_ebetac_list.dat");
	writeDataFile("generated_ebetac_list.dat", vdEbetacList);

	return vdEbetacList;
}

double getEbetac(std::string strFilename, InputData &sInput) 
{
	vector< vector<double> > vvdFesData = load2dDataFromFile(strFilename);
	double dNumerator = 0.0, dDenominator = 0.0;

	if (sInput.dGamma > 1e-6)  // Tempered
	{
		for (auto vdFesRow: vvdFesData)
		{
			double dExponent = - vdFesRow[sInput.nFreeEnergyColumn] / sInput.dkT;
			dNumerator += exp(dExponent);
			dDenominator += exp(dExponent / sInput.dGamma);
		}
	}
	else  // Not tempered
	{
		for (auto vdFesRow: vvdFesData)
			dNumerator += exp(- vdFesRow[sInput.nFreeEnergyColumn] / sInput.dkT);
		dDenominator = static_cast<double>(vvdFesData.size());
	}

	return dNumerator / dDenominator;
}

vector<double> loadEbetacList(InputData &sInput) 
{
	std::cout << "Loading the reversible work estimates... ";

	vector<double> vdEbetacList;
	if (!sInput.strFesDataFile.empty())
	{
		vdEbetacList = load1dDataFromFile(sInput.strFesDataFile);
		// Also change the FES files count to agree
		sInput.nFesFilesCount = vdEbetacList.size();
	}
	else  // Checked previously, so load from FES files is the only option
		vdEbetacList = loadFesFiles(sInput);

	std::cout << "Done" << std::endl;
	return vdEbetacList;
}

void assignInput(string strLine, InputData &sInput) 
{
	string strFirstInput, strSecondInput;

	size_t nPos = strLine.find_first_of('=');
	if (nPos == string::npos)
	{
		std::cerr << "The following line doesn't contain any assignment:" << std::endl;
		std::cerr << strLine << std::endl;
		exit (INVALID_INPUT);
	}
	else
	{
		strFirstInput = trimWhitespace(strLine.substr(0, nPos));
		strSecondInput = trimWhitespace(strLine.substr(nPos + 1));
	}

	// Assigning strings, simple
	if (!strFirstInput.compare("CoeffList"))			sInput.strCoeffList = strSecondInput;
	if (!strFirstInput.compare("ColvarFilePath"))		sInput.strColvarFilePath = strSecondInput;
	if (!strFirstInput.compare("FesDataFile"))			sInput.strFesDataFile = strSecondInput;
	if (!strFirstInput.compare("FesPrefix"))			sInput.strFesPrefix = strSecondInput;
	if (!strFirstInput.compare("LogFile"))				sInput.strLogFile = strSecondInput;
	if (!strFirstInput.compare("MetadColvarFile"))		sInput.strMetadColvarFile = strSecondInput;
	if (!strFirstInput.compare("NewColvarFile"))		sInput.strNewColvarFile = strSecondInput;
	if (!strFirstInput.compare("NewColvarFilePath"))	sInput.strNewColvarFilePath = strSecondInput;
	if (!strFirstInput.compare("PeriodicCoeffList"))	sInput.strPeriodicCoeffList = strSecondInput;

	// Assigning integers
	if (!strFirstInput.compare("AnnealingSteps"))		sInput.nAnnealingSteps = std::stoi(strSecondInput);
	if (!strFirstInput.compare("Barriers"))				sInput.nBarriers = std::stoi(strSecondInput);
	if (!strFirstInput.compare("FesFilesCount"))		sInput.nFesFilesCount = std::stoi(strSecondInput);
	if (!strFirstInput.compare("FreeEnergyColumn"))		sInput.nFreeEnergyColumn = std::stoi(strSecondInput);
	if (!strFirstInput.compare("Grid"))					sInput.nGrid = std::stoi(strSecondInput);	
	if (!strFirstInput.compare("LengthTarget"))			sInput.nLengthTarget = std::stoi(strSecondInput);
	if (!strFirstInput.compare("LengthTolerance"))		sInput.nLengthTolerance = std::stoi(strSecondInput);
	if (!strFirstInput.compare("LogEigenvalues"))		sInput.nLogEigenvalues = std::stoi(strSecondInput);
	if (!strFirstInput.compare("ScalingColumn"))		sInput.nScalingColumn = std::stoi(strSecondInput);
	if (!strFirstInput.compare("Seed"))					sInput.nSeed = std::stoi(strSecondInput);
	if (!strFirstInput.compare("TimeColumn"))			sInput.nTimeColumn = std::stoi(strSecondInput);

	// Assigning bools (via integers)
	if (!strFirstInput.compare("AlternativeInfinity"))	sInput.bAltInfinity = std::stoi(strSecondInput);
	if (!strFirstInput.compare("EuclideanDistance"))	sInput.bEuclidean = std::stoi(strSecondInput);
	if (!strFirstInput.compare("ConvertPeriodic"))		sInput.bConvertPeriodic = std::stoi(strSecondInput);
	if (!strFirstInput.compare("LimitIncorrect"))		sInput.bLimitIncorrect = std::stoi(strSecondInput);
	if (!strFirstInput.compare("LogPathEnergy"))		sInput.bLogEnergy = std::stoi(strSecondInput);
	if (!strFirstInput.compare("NoPath"))				sInput.bNoPath = std::stoi(strSecondInput);
	if (!strFirstInput.compare("Path2D"))				sInput.bPath2D = std::stoi(strSecondInput);
	if (!strFirstInput.compare("Random"))				sInput.bRandom = std::stoi(strSecondInput);
	if (!strFirstInput.compare("Rescale"))				sInput.bRescale = std::stoi(strSecondInput);
	if (!strFirstInput.compare("SmoothCount"))			sInput.bSmoothCount = std::stoi(strSecondInput);
	if (!strFirstInput.compare("SPath"))				sInput.bSPath = std::stoi(strSecondInput);
	if (!strFirstInput.compare("StrictOrder"))			sInput.bStrictOrder = std::stoi(strSecondInput);
	if (!strFirstInput.compare("TryAll"))				sInput.bTryAll = std::stoi(strSecondInput);

	// Assigning doubles
	if (!strFirstInput.compare("AnnealingkT"))			sInput.dAnnealingkT = std::stod(strSecondInput);
	if (!strFirstInput.compare("Cooling"))				sInput.dCooling = std::stod(strSecondInput);
	if (!strFirstInput.compare("Gamma"))				sInput.dGamma = std::stod(strSecondInput);
	if (!strFirstInput.compare("InfLimit"))				sInput.dInfLimit = std::stod(strSecondInput);
	if (!strFirstInput.compare("kT"))					sInput.dkT = std::stod(strSecondInput);
	if (!strFirstInput.compare("Threshold"))			sInput.dThreshold = std::stod(strSecondInput);
	if (!strFirstInput.compare("ZLimit"))				sInput.dZLimit = std::stod(strSecondInput);
	
	// Assigning vectors of integers
	if (!strFirstInput.compare("BiasColumns"))			sInput.vnBiasColumns = separateString(strSecondInput);
	if (!strFirstInput.compare("NewColumns"))			sInput.vnNewColumns = separateString(strSecondInput);
	if (!strFirstInput.compare("NewColumnsPath"))		sInput.vnNewColumnsPath = separateString(strSecondInput);
	if (!strFirstInput.compare("PeriodicColumns"))		sInput.vnPeriodicColumns = separateString(strSecondInput);
	if (!strFirstInput.compare("PeriodicColumnsPath"))	sInput.vnPeriodicColumnsPath = separateString(strSecondInput);

	// Assigning vectors of doubles
	if (!strFirstInput.compare("PeriodicRange"))		sInput.vdPeriodicRanges = separateStringDoubles(strSecondInput);
}

void checkInput(InputData &sInput) 
{
	bool bError = false;

	// Check input files
    if (sInput.strMetadColvarFile.empty())
    {
        std::cerr << "Both complete FES file and metadynamics Colvar file are unspecified, please specify one" << std::endl;
        bError = true;
    }
    if (sInput.strNewColvarFile.empty())
    {
        std::cerr << "Both complete FES file and new Colvar file are unspecified, please specify one" << std::endl;
        bError = true;
    }
    if (!sInput.strFesPrefix.empty() && !sInput.strFesDataFile.empty())
    {
        std::cerr << "Both FES prefix and FES data file have been specified, pick just one" << std::endl;
        bError = true;
    }
	if (sInput.strNewColvarFilePath.empty())
	{
		std::cerr << "New Colvar file for path is unspecified, please specify one" << std::endl;
		bError = true;
	}
	if (sInput.strFesPrefix.empty() && sInput.strFesDataFile.empty())
	{
		std::cerr << "Both FES data file and FES prefix are unspecified, please specify one of them" << std::endl;
		bError = true;
	}
	if (!sInput.strCoeffList.empty() && sInput.bConvertPeriodic && sInput.strPeriodicCoeffList.empty())
	{
		std::cerr << "Coefficient list for periodic variables is unspecified but needed, please include it" << std::endl;
		bError = true;
	}

	// Check vectors
	if (sInput.vnBiasColumns.empty())
	{
		std::cerr << "No bias columns have been specified, please specify at least one" << std::endl;
		bError = true;
	}
	if (sInput.vnNewColumns.empty())
	{
		std::cerr << "No new CV columns have been specified, please specify at least one" << std::endl;
		bError = true;
	}
	if (sInput.vnNewColumnsPath.empty())
	{
		std::cerr << "No new CV columns for path have been specified, please specify at least one" << std::endl;
		bError = true;
	}
	if (sInput.vdPeriodicRanges.size() != sInput.vnPeriodicColumns.size())
	{
		std::cerr << "The number of periodic variables must match the number of ranges provided" << std::endl;
		bError = true;
	}
	if (sInput.vnPeriodicColumnsPath.size() != sInput.vnPeriodicColumns.size())
	{
		std::cerr << "The number of periodic variables for the path and metadynamics must be the same" << std::endl;
		bError = true;
	}
	if (sInput.vnNewColumnsPath.size() != sInput.vnNewColumns.size())
	{
		std::cerr << "The number of collective variables for the path and metadynamics must be the same" << std::endl;
		bError = true;
	}

	// Check stuff only relevant for FES file loading
	if (!sInput.strFesPrefix.empty())
	{
		if (sInput.nFesFilesCount == -1)
		{
			std::cerr << "Number of FES files is not specified or invalid, please specify it" << std::endl;
			bError = true;
		}
		if (sInput.nFreeEnergyColumn == -1)
		{
			std::cerr << "The free energy column is not specified or invalid, please specify it" << std::endl;
			bError = true;
		}
		if (sInput.dGamma < 0.0)
		{
			std::cerr << "The biasfactor is not specified or invalid, please specify it" << std::endl;
			bError = true;
		}
	}

	// Check all the other integers
	if (sInput.nGrid == -1)
	{
		std::cerr << "The number of bins for reweighting is not specified or invalid, please specify it" << std::endl;
		bError = true;
	}
	if (sInput.nLengthTarget == -1)
	{
		std::cerr << "The number of snapshots is not specified or invalid, please specify it" << std::endl;
		bError = true;
	}
	if (sInput.nAnnealingSteps == -1)
	{
		std::cerr << "The number of annealing steps is not specified or invalid, please specify it" << std::endl;
		bError = true;
	}

	// Check all the other doubles
	if (sInput.dkT < 0.0)
	{
		std::cerr << "The value of kT is not specified or invalid, please specify it" << std::endl;
		bError = true;
	}
	if (sInput.dAnnealingkT < 0.0)
	{
		std::cerr << "The value of initial kT for annealing is not specified or invalid, please specify it" << std::endl;
		bError = true;
	}

	if (bError)
		exit (INVALID_INPUT);
}

InputData loadInput(std::string strInputFile) 
{
	InputData sInput;

	std::ifstream InputStream(strInputFile.c_str());

	if (!InputStream)
	{
		std::cerr << "Couldn't open file " << strInputFile << std::endl;
		exit (FILE_NOT_FOUND);
	}

	while (InputStream)
	{
		std::string strLine;
		getline(InputStream, strLine);

		// Ignore comment or empty lines
		if (!std::strncmp(strLine.c_str(), "#", 1) || isEmpty(strLine)) 
			continue;

		// Remove everything after #
		size_t nPos = strLine.find_first_of('#');
		if (nPos != std::string::npos)
			strLine.erase(nPos);

		assignInput(strLine, sInput);
	}

	return sInput;
}

vector<double> load1dDataFromFile(std::string strFilename) 
{
	vector<double> vdData;

	std::ifstream InputStream(strFilename.c_str());

	if (!InputStream)
	{
		std::cerr << "Couldn't open file " << strFilename << std::endl;
		exit (FILE_NOT_FOUND);
	}

	while (InputStream)
	{
		std::string strLine;
		getline(InputStream, strLine);
		// Ignore comment or empty lines
		if (!std::strncmp(strLine.c_str(), "#", 1) || isEmpty(strLine)) 
			continue;

		double dValue = std::atof(strLine.c_str());
		vdData.push_back(dValue);
	}

	return vdData;
}

vector< vector<double> > load2dDataFromFile(std::string strFilename) 
{
	vector< vector<double> > vvdData;

	std::ifstream InputStream(strFilename);

	if (!InputStream)
	{
		std::cerr << "Couldn't open file " << strFilename << std::endl;
		exit (FILE_NOT_FOUND);
	}

	while (InputStream)
	{
		std::string strLine;
		getline(InputStream, strLine);
		// Ignore comment or empty lines
		if (!std::strncmp(strLine.c_str(), "#", 1) || isEmpty(strLine)) 
			continue;

		std::stringstream StringStream;
		StringStream << strLine;

		vector<double> vdDataLine;
		while (!StringStream.eof())
		{
			std::string strValue;
			StringStream >> strValue;

			if (!strValue.empty())
			{
				double dValue = std::atof(strValue.c_str());
				vdDataLine.push_back(dValue);
			}
		}
		vvdData.push_back(vdDataLine);
	}

	return vvdData;
}
