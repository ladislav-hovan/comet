/*
 * testing.cpp
 *
 *  Created on: 13 Apr 2016
 *      Author: ladislav
 */

#include "testing.h"

void addRandomColumn(vector< vector<double> > &vvdData) {
	std::random_device Random;  // Random seed
	std::mt19937 Mersenne(Random());
	std::uniform_real_distribution<double> Distribution(0.0, 1.0);

	for (auto &vdRow: vvdData)
		vdRow.push_back(Distribution(Mersenne));
}

void addConstantColumn(vector< vector<double> > &vvdData, double dConstant) {
	for (auto &vdRow: vvdData)
		vdRow.push_back(dConstant);
}


