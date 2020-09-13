// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT
//
// Derived from https://github.com/osqzss/gps-sdr-sim (MIT Licence):
// Copyright (c) 2015-2020 Takuji Ebinuma

#pragma once

#include <Eigen/Core>

#include <cmath>

/*! \brief Subtract two vectors of double
 *  \param[out] y Result of subtraction
 *  \param[in] x1 Minuend of subtracion
 *  \param[in] x2 Subtrahend of subtracion
 */
static inline void subVect(double *y, const double *x1, const double *x2)
{
	y[0] = x1[0]-x2[0];
	y[1] = x1[1]-x2[1];
	y[2] = x1[2]-x2[2];

	return;
}

/*! \brief Subtract two vectors of double
 *  \param[out] y Result of subtraction
 *  \param[in] x1 Minuend of subtracion
 *  \param[in] x2 Subtrahend of subtracion
 */
static inline void subVect(double *y, const double *x1, const Eigen::Vector3d &x2)
{
	y[0] = x1[0]-x2.coeff(0);
	y[1] = x1[1]-x2.coeff(1);
	y[2] = x1[2]-x2.coeff(2);
}

/*! \brief Compute Norm of Vector
 *  \param[in] x Input vector
 *  \returns Length (Norm) of the input vector
 */
static inline double normVect(const double *x)
{
	return std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
}

/*! \brief Compute dot-product of two vectors
 *  \param[in] x1 First multiplicand
 *  \param[in] x2 Second multiplicand
 *  \returns Dot-product of both multiplicands
 */
static inline double dotProd(const double *x1, const double *x2)
{
	return(x1[0]*x2[0]+x1[1]*x2[1]+x1[2]*x2[2]);
}
