// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT
//
// Derived from https://github.com/osqzss/gps-sdr-sim (MIT Licence):
// Copyright (c) 2015-2020 Takuji Ebinuma

#pragma once

#include <array>

class AntennaGain
{
  public:
      AntennaGain();

      double ComputeAntennaGain(double azel_rad[2]) const;

  private:
      // Receiver antenna attenuation for boresight angle = 0:5:180 [deg] as
      // natural units. The horizon is 0deg and 180deg. Zenith is 90deg.
      std::array<double, 37> gain_pattern_boresight_;
};
