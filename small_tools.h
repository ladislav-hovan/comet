/*
 * small_tools.h
 *
 *  Created on: 10 Mar 2016
 *      Author: ladislav
 */

#ifndef SMALL_TOOLS_H_
#define SMALL_TOOLS_H_

#include <string>
#include <cstring>
#include <cctype>
#include <iostream>
#include <cmath>
#include "input_data.h"
#include "error_codes.h"

using std::string;

bool isEmpty(string strInput);
void handleSeparated(vector<int> &vnSequence, bool &bRange, string &strMember, int nRangeStart);
vector<int> separateString(string strSequence);
vector<double> separateStringDoubles(string strSequence);
string trimWhitespace(string strInput);
double setDecimal(double dValue, int nFigures);
double setSignificant(double dValue, int nFigures);
void setSignificant(vector<double>& vdValues, int nFigures);

#endif /* SMALL_TOOLS_H_ */
