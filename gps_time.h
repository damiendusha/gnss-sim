// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT
//
// Derived from https://github.com/osqzss/gps-sdr-sim (MIT Licence):
// Copyright (c) 2015-2020 Takuji Ebinuma

#pragma once

#include "gpssim.h"

void date2gps(const datetime_t *t, gpstime_t *g);
void gps2date(const gpstime_t *g, datetime_t *t);
gpstime_t incGpsTime(gpstime_t g0, double dt);
double subGpsTime(gpstime_t g1, gpstime_t g0);
