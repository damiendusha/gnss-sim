// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT

#include <random>

class NoiseGenerator
{
  public:
    NoiseGenerator(double signal_gain, double noise_gain, double output_gain)
        : signal_gain_(signal_gain)
        , output_gain_(output_gain)
        , awgn_(kZeroMean, noise_gain)
    {
    }
    
    double ScaleAndAddNoise(double sample)
    {
        const double noise = awgn_(generator_);
        return (signal_gain_ * sample + noise) * output_gain_;
    }
      
  private:
    static constexpr double kZeroMean = 0.0;

    const double signal_gain_;
    const double output_gain_;

    std::default_random_engine generator_;
    std::normal_distribution<double> awgn_;
};
