/*
 * small_tools.cpp
 *
 *  Created on: 11 Mar 2016
 *      Author: ladislav
 */

#include "small_tools.h"

bool isEmpty(string strInput) 
{
	int nCount = 0;

	while (strInput[nCount])
	{
		if (!std::isspace(strInput[nCount])) return false;
		++nCount;
	}

	return true;
}

void handleSeparated(vector<int> &vnSequence, bool &bRange, string &strMember, int nRangeStart)
{
	int nMember = std::stoi(strMember);

	if (bRange)
	{
		if (nMember < nRangeStart)
		{
			std::cerr << "The range in the input file is invalid: " << nRangeStart << "-" << nMember << std::endl;
			exit(INVALID_INPUT);
		}
		else
		{
			while (nRangeStart <= nMember)
				vnSequence.push_back(nRangeStart++);
		}

		bRange = false;
	}
	else
		vnSequence.push_back(nMember);

	strMember.clear();
}

vector<int> separateString(string strSequence)
{
	vector<int> vnSequence;
	string strMember;
	int nCount = 0, nRangeStart = 0;
	bool bRange = false;

	while (strSequence[nCount])
	{
		if (!strncmp(&strSequence[nCount], ",", 1))
			handleSeparated(vnSequence, bRange, strMember, nRangeStart);
		else if (!strncmp(&strSequence[nCount], "-", 1))
		{
			nRangeStart = std::stoi(strMember);
			bRange = true;
			strMember.clear();
		}
		else
			strMember += strSequence[nCount];

		++nCount;
	}

	handleSeparated(vnSequence, bRange, strMember, nRangeStart);

	return vnSequence;
}

vector<double> separateStringDoubles(string strSequence)
{
	vector<double> vdSequence;
	string strMember;
	int nCount = 0;

	while (strSequence[nCount])
	{
		if (!strncmp(&strSequence[nCount], ",", 1))
		{
			double dMember = std::stod(strMember);
			vdSequence.push_back(dMember);
			strMember.clear();
		}
		else
			strMember += strSequence[nCount];

		++nCount;
	}

	double dMember = std::stod(strMember);
	vdSequence.push_back(dMember);
	strMember.clear();

	return vdSequence;
}

string trimWhitespace(string strInput)
{
	size_t nFirst = strInput.find_first_not_of(" \t");
	size_t nSecond = strInput.find_last_not_of(" \t");

	string strOutput = strInput.substr(nFirst, nSecond - nFirst + 1);

	return strOutput;
}

double setDecimal(double dValue, int nFigures)
{
	double dShift = std::pow(10, nFigures);

	return std::round(dValue * dShift) / dShift;
}

double setSignificant(double dValue, int nFigures)
{
	int nLog = std::ceil(std::log10(dValue));
	double dShift = std::pow(10, nFigures - nLog);

	return std::round(dValue * dShift) / dShift;
}

void setSignificant(vector<double> &vdValues, int nFigures)
{
	for (double& dValue : vdValues)
	{
		int nLog = std::ceil(std::log10(dValue));
		double dShift = std::pow(10, nFigures - nLog);

		dValue = std::round(dValue * dShift) / dShift;
	}
}
