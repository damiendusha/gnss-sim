// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT
//
// Derived from https://github.com/osqzss/gps-sdr-sim (MIT Licence):
// Copyright (c) 2015-2020 Takuji Ebinuma

#include "ionosphere.h"

namespace {

static constexpr double SECONDS_IN_DAY = 86400.0;
static constexpr double SPEED_OF_LIGHT = 2.99792458e8;
static constexpr double PI = 3.1415926535898; // GPS-ICD-200 !!!
}

double ionosphericDelay(const ionoutc_t &ionoutc, gpstime_t g, double *llh, 
                        const AzimuthElevation &azel)
{
	if (ionoutc.enable == false)
		return 0; // No ionospheric delay

	const double E = azel.elevation_rad() / PI;
	const double phi_u = llh[0]/PI;
	const double lam_u = llh[1]/PI;

	// Obliquity factor
	const double F = 1.0 + 16.0*pow((0.53 - E),3.0);

    double iono_delay = 0.0;
	if (!ionoutc.valid) {
		iono_delay = F*5.0e-9*SPEED_OF_LIGHT;
    }
	else
	{
		// Earth's central angle between the user position and the earth projection of
		// ionospheric intersection point (semi-circles)
		const double psi = 0.0137/(E + 0.11) - 0.022;
		
		// Geodetic latitude of the earth projection of the ionospheric intersection point
		// (semi-circles)
		double phi_i = phi_u + psi * std::cos(azel.azimuth_rad());
		if(phi_i>0.416)
			phi_i = 0.416;
		else if(phi_i<-0.416)
			phi_i = -0.416;

		// Geodetic longitude of the earth projection of the ionospheric intersection point
		// (semi-circles)
		const double lam_i = lam_u + psi* std::sin(azel.azimuth_rad()) / std::cos(phi_i*PI);

		// Geomagnetic latitude of the earth projection of the ionospheric intersection
		// point (mean ionospheric height assumed 350 km) (semi-circles)
		const double phi_m = phi_i + 0.064 * std::cos((lam_i - 1.617)*PI);
		const double phi_m2 = phi_m*phi_m;
		const double phi_m3 = phi_m2*phi_m;

		double AMP = ionoutc.alpha0 + ionoutc.alpha1*phi_m
			+ ionoutc.alpha2*phi_m2 + ionoutc.alpha3*phi_m3;
		if (AMP < 0.0)
			AMP = 0.0;

		double PER = ionoutc.beta0 + ionoutc.beta1*phi_m
			+ ionoutc.beta2*phi_m2 + ionoutc.beta3*phi_m3;
		if (PER < 72000.0)
			PER = 72000.0;

		// Local time (sec)
		double t = SECONDS_IN_DAY/2.0*lam_i + g.sec;
		while(t>=SECONDS_IN_DAY)
			t -= SECONDS_IN_DAY;
		while(t<0)
			t += SECONDS_IN_DAY;

		// Phase (radians)
		const double X = 2.0*PI*(t - 50400.0)/PER;

		if (std::fabs(X)<1.57)
		{
			const double X2 = X*X;
			const double X4 = X2*X2;
			iono_delay = F*(5.0e-9 + AMP*(1.0 - X2/2.0 + X4/24.0))*SPEED_OF_LIGHT;
		}
		else {
			iono_delay = F*5.0e-9*SPEED_OF_LIGHT;
        }
	}

	return iono_delay;
}

