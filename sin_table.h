// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT

#pragma once

#include <functional>
#include <vector>
#include <iostream>
#include <cmath>

template <typename T, int kTableLength>
class PeriodicFunctionTable
{
  public:
    explicit PeriodicFunctionTable(const std::function<double(double)> &func, 
                                   double amplitude)
    {
        const double scale = 2.0 * M_PI / kTableLength;
        for (int i = 0; i < kTableLength; ++i) {
            const double angle_rad = scale * i;
            table_[i] = amplitude * func(angle_rad);
        }
    }

    // Normalised phase is the phase in circles. The mapping is as follows:
    // 0    =   0 deg
    // 0.25 =  90 deg
    // 0.5  = 180 deg 
    // 0.75 = 270 deg
    // 1.0  = 360 deg
    double LookupValue(double normalised_phase) const
    {
        const double raw_index = normalised_phase * table_.size();
        const int index = std::floor(raw_index);   
        return table_[index];
    }

    void PrintTable() const 
    {
        for (int i = 0 ; i < table_.size(); ++i) {
            std::cout << "[" << i << "] = " << table_[i] << std::endl;
        }
    }

  private:
    std::array<T, kTableLength> table_;
};
