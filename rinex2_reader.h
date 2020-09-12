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
#include "gps_ephem.h"

#include <string>

/*! \brief Read Ephemeris data from the RINEX Navigation file */
/*  \param[out] eph Array of Output SV ephemeris data
 *  \param[in] fname File name of the RINEX file
 *  \returns Number of sets of ephemerides in the file
 */
int readRinexNavAll(ephem_t eph[][MAX_SAT], ionoutc_t *ionoutc, 
                    const std::string &rinex_filename);

