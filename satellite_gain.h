// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT

#pragma once

#include "antenna_gain.h"

#include <memory>
#include <array>

class SatelliteGain
{
  public:
    explicit SatelliteGain(double scale = 1.0) : scale_(scale) {}
    virtual ~SatelliteGain() {}
      
    virtual double ComputeGain(double range_m, double azel_rad[2])
    {
        return scale_ * InternalComputeGain(range_m, azel_rad);
    }

  protected:
    const double scale_;

    virtual double InternalComputeGain(double range_m, double azel_rad[2]) const = 0;
};

class SatelliteGainRangeRxAntenna : public SatelliteGain
{
  public:

  protected:
    double InternalComputeGain(double range_m, double azel_rad[2]) const override 
    {
		const double path_loss = 20200000.0 / range_m;
		const double ant_gain = antenna_gain_.ComputeAntennaGain(azel_rad);

		// Signal gain
		return path_loss * ant_gain;
    }

  private:
    AntennaGain antenna_gain_;
};

class SatelliteGainConstant : public SatelliteGain
{
  public:
    explicit SatelliteGainConstant(const double gain) : gain_(gain) {}
      
  protected:
    double InternalComputeGain(double range_m, double azel_rad[2]) const override 
    {
        return gain_;
    }

  private:
    const double gain_;
};

class ConstellationGain
{
  public:
    ConstellationGain();

    double ComputeGain(int prn, double range_m, double azel_rad[2]) const;

    void SetSatelliteToConstantGain(int prn, double gain);
    
  private:
    static int IndexFromPrn(int prn) { return prn - 1; }

    std::array<std::unique_ptr<SatelliteGain>, 32> prn_to_gain_;
};
