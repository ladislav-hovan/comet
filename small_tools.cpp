/*
 * small_tools.cpp
 *
 *  Created on: 11 Mar 2016
 *      Author: ladislav
 */

#include "small_tools.h"

bool isEmpty(std::string strInput) 
{
	int nCount = 0;

	while (strInput[nCount])
	{
		if (!std::isspace(strInput[nCount])) return false;
		++nCount;
	}

	return true;
}

void handleSeparated(vector<int> &vnSequence, bool &bRange, std::string &strMember, int nRangeStart)
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

vector<int> separateString(std::string strSequence)
{
	vector<int> vnSequence;
	std::string strMember;
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

vector<double> separateStringDoubles(std::string strSequence)
{
	vector<double> vdSequence;
	std::string strMember;
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

std::string trimWhitespace(std::string strInput)
{
	size_t nFirst = strInput.find_first_not_of(" \t");
	size_t nSecond = strInput.find_last_not_of(" \t");

	std::string strOutput = strInput.substr(nFirst, nSecond - nFirst + 1);

	return strOutput;
}
