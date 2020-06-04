// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT

#include "satellite_gain.h"

double ConstellationGain::ComputeGain(int prn, double range_m, double azel_rad[2]) const
{
    const int index = IndexFromPrn(prn);
    return prn_to_gain_[index]->ComputeGain(range_m, azel_rad);
}

void ConstellationGain::SetSatelliteToConstantGain(int prn, double gain)
{
    int index = IndexFromPrn(prn);
    prn_to_gain_[index] = std::make_unique<SatelliteGainConstant>(gain);
}

ConstellationGain::ConstellationGain()
{
    for (int i = 0; i < prn_to_gain_.size(); ++i)
    {
        prn_to_gain_[i] = std::make_unique<SatelliteGainRangeRxAntenna>();
    }
}
