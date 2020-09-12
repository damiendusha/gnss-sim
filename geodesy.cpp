// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT
//
// Derived from https://github.com/osqzss/gps-sdr-sim (MIT Licence):
// Copyright (c) 2015-2020 Takuji Ebinuma

#include "geodesy.h"
#include "gpssim.h"
#include "gps_math.h"

#include <cmath>

GeodeticPosition xyz2llh(const double *xyz)
{
	double rho2,dz,zdz,nh,slat,n,dz_new;

	constexpr double a = WGS84_RADIUS;
	constexpr double e = WGS84_ECCENTRICITY;
	constexpr double e2 = e*e;

    constexpr double eps = 1.0e-3;
	if (normVect(xyz) < eps)
	{
        // Invalid ECEF vector
        return GeodeticPosition::FromRadians(0.0, 0.0, -a);
	}

	const double x = xyz[0];
	const double y = xyz[1];
	const double z = xyz[2];

	rho2 = x*x + y*y;
	dz = e2*z;

	while (1)
	{
		zdz = z + dz;
		nh = std::sqrt(rho2 + zdz*zdz);
		slat = zdz / nh;
		n = a / std::sqrt(1.0-e2*slat*slat);
		dz_new = n*e2*slat;

		if (fabs(dz-dz_new) < eps)
			break;

		dz = dz_new;
	}

    const double lat_rad = std::atan2(zdz, std::sqrt(rho2));
    const double lon_rad = std::atan2(y, x);
    const double height_m = nh - n;

    return GeodeticPosition::FromRadians(lat_rad, lon_rad, height_m);
}

void llh2xyz(const GeodeticPosition &llh, double *xyz)
{
	constexpr double a = WGS84_RADIUS;
	constexpr double e = WGS84_ECCENTRICITY;
	constexpr double e2 = e*e;

	const double clat = std::cos(llh.latitude_rad());
	const double slat = std::sin(llh.latitude_rad());
	const double clon = std::cos(llh.longitude_rad());
	const double slon = std::sin(llh.longitude_rad());
	const double d = e*slat;

	const double n = a / std::sqrt(1.0-d*d);
	const double nph = n + llh.height_m();

	const double tmp = nph*clat;
	xyz[0] = tmp*clon;
	xyz[1] = tmp*slon;
	xyz[2] = ((1.0-e2)*n + llh.height_m())*slat;

	return;
}

void ltcmat(const GeodeticPosition &llh, double t[3][3])
{
	const double slat = std::sin(llh.latitude_rad());
	const double clat = std::cos(llh.latitude_rad());
	const double slon = std::sin(llh.longitude_rad());
	const double clon = std::cos(llh.longitude_rad());

	t[0][0] = -slat*clon;
	t[0][1] = -slat*slon;
	t[0][2] = clat;
	t[1][0] = -slon;
	t[1][1] = clon;
	t[1][2] = 0.0;
	t[2][0] = clat*clon;
	t[2][1] = clat*slon;
	t[2][2] = slat;
}

void ecef2neu(const double *xyz, double t[3][3], double *neu)
{
	neu[0] = t[0][0]*xyz[0] + t[0][1]*xyz[1] + t[0][2]*xyz[2];
	neu[1] = t[1][0]*xyz[0] + t[1][1]*xyz[1] + t[1][2]*xyz[2];
	neu[2] = t[2][0]*xyz[0] + t[2][1]*xyz[1] + t[2][2]*xyz[2];

	return;
}

AzimuthElevation neu2azel(const double *neu)
{
	double azimuth_rad = std::atan2(neu[1],neu[0]);
	if (azimuth_rad < 0.0)
		azimuth_rad += 2.0*PI;
    
	const double ne = std::hypot(neu[0], neu[1]);
	const double elevation_rad = std::atan2(neu[2], ne);

    return AzimuthElevation::FromRadians(azimuth_rad, elevation_rad);
}
