// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT
//
// Derived from https://github.com/osqzss/gps-sdr-sim (MIT Licence):
// Copyright (c) 2015-2020 Takuji Ebinuma

#ifndef GPS_CHANNEL_H_
#define GPS_CHANNEL_H_

#include "gpssim.h"
#include "gps_ephem.h"
#include "ionosphere.h"

#include <array>

/*! \brief Structure representing a Channel */
class GpsChannel
{
  public:
    GpsChannel() = default;
    GpsChannel(int in_prn);

    /// \brief PRN for the satellite. Zero if the channel is not allocated.
    int prn = 0;

    // Current carrier phase in cycles (i.e. value of 0 to 1).
    double carrier_phase_cycles;

    /// \brief GPS time at the start of the full data frame (i.e. 30s boundary).
    gpstime_t dataframe_reference_time;

	AzimuthElevation azel;
	range_t rho0;

    /// \brief Returns true if the channel is enabled.
    bool IsEnabled() const { return prn > 0; }
    
    void ComputeCodePhase(const range_t &rho1, const double dt);

    double current_code_data_symbol() const { 
        return current_code_data_symbol_;
    }

    void UpdateCarrierPhase(const double delta_t_s) 
    {
        carrier_phase_cycles += f_carr * delta_t_s;
        if (carrier_phase_cycles >= 1.0)
            carrier_phase_cycles -= 1.0;
        else if (carrier_phase_cycles < 0.0)
            carrier_phase_cycles += 1.0;
    }

    void UpdateCodePhase(const double delta_t_s)
    {
        code_phase += f_code * delta_t_s;

        if (code_phase >= CA_SEQ_LEN) {
            code_phase -= CA_SEQ_LEN;

            initial_code++;
            if (initial_code >= 20) // 20 C/A codes = 1 navigation data bit
            {
                initial_code = 0;
                
                initial_bit++;
                if (initial_bit >= 30) // 30 navigation data bits = 1 word
                {
                    initial_bit = 0;
                    initial_word++;
                }

                // Set new navigation data bit
                UpdateDataBit();
            }
        }
        
        UpdateCodeChip();
    }
    
    void GenerateNavMsg(gpstime_t g, int init);

    void SetEphemeris(const ephem_t &eph, const ionoutc_t &ionoutc);


  private:
    /// \brief C/A code sequence. All values are either -1 or 1.
    std::array<int, CA_SEQ_LEN> code_sequence_;
    
    int current_data_bit_ = 0;
    int current_code_chip_ = 0;
    double current_code_data_symbol_ = 0;
    
    int initial_word;	/*!< initial word */
    int initial_bit;	/*!< initial bit */
    int initial_code;	/*!< initial code */
    
    /// \brief Carrier frequency.
    double f_carr;

    /// \brief Code frequency.
    double f_code;
    
    double code_phase; /*< Code phase */

    /// \brief Current subframe.
    unsigned long sbf[5][N_DWRD_SBF];

    /// \brief Data words of the subframe.
    unsigned long dwrd[N_DWRD];
    
    int ComputeDataBit() const {
        return (int)((dwrd[initial_word] >> (29 - initial_bit)) & 0x1UL)*2 - 1;
    }

    int ComputeCodeChip() const {
        return code_sequence_[(int) code_phase];
    }
    
    void UpdateCodeChip() {
        current_code_chip_ = ComputeCodeChip();
        current_code_data_symbol_ = current_code_chip_ * current_data_bit_;
    }

    void UpdateDataBit() {
        current_data_bit_ = ComputeDataBit();
        current_code_data_symbol_ = current_code_chip_ * current_data_bit_;
    }
};

#endif // GPS_CHANNEL_H_

