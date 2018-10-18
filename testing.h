/*
 * testing.h
 *
 *  Created on: 13 Apr 2016
 *      Author: ladislav
 */

#include <vector>
#include <random>

using std::vector;

#ifndef TESTING_H_
#define TESTING_H_

void addRandomColumn(vector< vector<double> > &vvdData);
void addConstantColumn(vector< vector<double> > &vvdData, double dConstant);

#endif /* TESTING_H_ */
