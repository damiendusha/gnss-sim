// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT
//
// Derived from https://github.com/osqzss/gps-sdr-sim (MIT Licence):
// Copyright (c) 2015-2020 Takuji Ebinuma

#include "gps_time.h"

#include <cmath>
#include <array>

namespace {

static constexpr double SECONDS_IN_WEEK = 604800.0;
static constexpr double SECONDS_IN_HALF_WEEK = 302400.0;
static constexpr double SECONDS_IN_DAY = 86400.0;
static constexpr double SECONDS_IN_HOUR = 3600.0;
static constexpr double SECONDS_IN_MINUTE = 60.0;

}

/*! \brief Convert a UTC date into a GPS date
 *  \param[in] t input date in UTC form
 *  \param[out] g output date in GPS form
 */
gpstime_t date2gps(const datetime_t &t)
{
	constexpr std::array<int, 12> doy = {
        0,31,59,90,120,151,181,212,243,273,304,334};

	const int ye = t.y - 1980;

	// Compute the number of leap days since Jan 5/Jan 6, 1980.
	int lpdays = ye/4 + 1;
	if ((ye%4)==0 && t.m<=2)
		lpdays--;

	// Compute the number of days elapsed since Jan 5/Jan 6, 1980.
	const int de = ye*365 + doy[t.m-1] + t.d + lpdays - 6;

	// Convert time to GPS weeks and seconds.
    gpstime_t g;
	g.week = de / 7;
	g.sec = (double)(de%7)*SECONDS_IN_DAY + t.hh*SECONDS_IN_HOUR 
		+ t.mm*SECONDS_IN_MINUTE + t.sec;

	return g;
}

datetime_t gps2date(const gpstime_t &g)
{
	// Convert Julian day number to calendar date
	const int c = (int)(7*g.week + std::floor(g.sec/86400.0)+2444245.0) + 1537;
	const int d = (int)((c-122.1)/365.25);
	const int e = 365*d + d/4;
	const int f = (int)((c-e)/30.6001);

    datetime_t t;
	t.d = c - e - (int)(30.6001*f);
	t.m = f - 1 - 12*(f/14);
	t.y = d - 4715 - ((7 + t.m)/10);

	t.hh = ((int)(g.sec/3600.0))%24;
	t.mm = ((int)(g.sec/60.0))%60;
	t.sec = g.sec - 60.0*floor(g.sec/60.0);

	return t;
}

double subGpsTime(gpstime_t g1, gpstime_t g0)
{
	double dt;

	dt = g1.sec - g0.sec;
	dt += (double)(g1.week - g0.week) * SECONDS_IN_WEEK;

	return(dt);
}

gpstime_t incGpsTime(gpstime_t g0, double dt)
{
	gpstime_t g1;

	g1.week = g0.week;
	g1.sec = g0.sec + dt;

	g1.sec = round(g1.sec*1000.0)/1000.0; // Avoid rounding error

	while (g1.sec>=SECONDS_IN_WEEK)
	{
		g1.sec -= SECONDS_IN_WEEK;
		g1.week++;
	}

	while (g1.sec<0.0)
	{
		g1.sec += SECONDS_IN_WEEK;
		g1.week--;
	}

	return(g1);
}
