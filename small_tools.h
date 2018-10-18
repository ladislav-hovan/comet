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
#include "input_data.h"
#include "error_codes.h"

bool isEmpty(std::string strInput);
void handleSeparated(vector<int> &vnSequence, bool &bRange, std::string &strMember, int nRangeStart);
vector<int> separateString(std::string strSequence);
vector<double> separateStringDoubles(std::string strSequence);
std::string trimWhitespace(std::string strInput);

#endif /* SMALL_TOOLS_H_ */
