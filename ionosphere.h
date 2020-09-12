// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT
//
// Derived from https://github.com/osqzss/gps-sdr-sim (MIT Licence):
// Copyright (c) 2015-2020 Takuji Ebinuma

#ifndef IONOSPHERE_H_
#define IONOSPHERE_H_

#include "geodesy.h"
#include "gps_time.h"

struct ionoutc_t
{
	bool enable;
	bool valid;
	double alpha0,alpha1,alpha2,alpha3;
	double beta0,beta1,beta2,beta3;
	double A0,A1;
	int dtls,tot,wnt;
	int dtlsf,dn,wnlsf;
};

double ionosphericDelay(const ionoutc_t &ionoutc, gpstime_t g, double *llh, 
                        const AzimuthElevation &azel);


#endif  // IONOSPHERE_H_
