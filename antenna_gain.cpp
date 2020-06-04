// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT
//
// Derived from https://github.com/osqzss/gps-sdr-sim (MIT Licence):
// Copyright (c) 2015-2020 Takuji Ebinuma

#include "antenna_gain.h"
#include "gpssim.h"

#include <cmath>

namespace
{
    // Receiver antenna attenuation in dB for boresight angle = 0:5:180 [deg]
    constexpr std::array<double, 37> kDefaultAntennaGain_db = {
         0.00,  0.00,  0.22,  0.44,  0.67,  1.11,  1.56,  2.00,  2.44,  2.89,  3.56,  4.22,
         4.89,  5.56,  6.22,  6.89,  7.56,  8.22,  8.89,  9.78, 10.67, 11.56, 12.44, 13.33,
        14.44, 15.56, 16.67, 17.78, 18.89, 20.00, 21.33, 22.67, 24.00, 25.56, 27.33, 29.33,
        31.56
    };
}   // namespace

AntennaGain::AntennaGain()
{
    for (int i = 0; i < 37; i++) 
    {
        gain_pattern_boresight_[i] = std::pow(10.0, -kDefaultAntennaGain_db[i] / 20.0);
    }
}

double AntennaGain::ComputeAntennaGain(double azel_rad[2]) const
{
    // Receiver antenna gain
	const int ibs = (int)((90.0-azel_rad[1]*R2D)/5.0); // convert elevation to boresight.
    return gain_pattern_boresight_[ibs];
}
