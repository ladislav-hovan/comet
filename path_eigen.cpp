/*
 * path_eigen.cpp
 *
 *  Created on: 29 Sep 2016
 *      Author: ladislav
 */

// This function is separate in hopes that it won't need to be recompiled that often

#include "path.h"

// Not included in the header, just used for sorting
bool descending(double dFirst, double dSecond) 
{
	return dFirst > dSecond;
}

void Path::calculateSpectralGap() 
{
	//logPathEnergies();
	populatekMatrix();
	Eigen::EigenSolver<Eigen::MatrixXd> Solver(m_kMatrix, false);  // False since we do not need eigenvectors
	m_vdEigenvalues.clear();
	for (unsigned int nCount = 0; nCount < m_vnSnapshots.size(); ++nCount)
		// Use only real part, imaginary should be zero but not checking
		m_vdEigenvalues.push_back(Solver.eigenvalues()[nCount].real());  
	std::sort(m_vdEigenvalues.begin(), m_vdEigenvalues.end(), descending);
	//writeDataFile("eigenvalues.dat", m_vdEigenvalues);
	determineSpectralGap();
}
