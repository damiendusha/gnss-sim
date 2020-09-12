// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT
//
// Derived from https://github.com/osqzss/gps-sdr-sim (MIT Licence):
// Copyright (c) 2015-2020 Takuji Ebinuma

#pragma once

/*! \brief Structure representing GPS time */
struct gpstime_t
{
    int week;	/*!< GPS week number (since January 1980) */
    double sec; 	/*!< second inside the GPS \a week */
};

/*! \brief Structure representing UTC time */
struct datetime_t
{
	int y; 		/*!< Calendar year */
	int m;		/*!< Calendar month */
	int d;		/*!< Calendar day */
	int hh;		/*!< Calendar hour */
	int mm;		/*!< Calendar minutes */
	double sec;	/*!< Calendar seconds */
};

gpstime_t date2gps(const datetime_t &t);
datetime_t gps2date(const gpstime_t &g);
gpstime_t incGpsTime(gpstime_t g0, double dt);
double subGpsTime(gpstime_t g1, gpstime_t g0);
