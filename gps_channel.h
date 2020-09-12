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

#include <array>

/*! \brief Structure representing a Channel */
class GpsChannel
{
  public:
    GpsChannel() = default;
    GpsChannel(int in_prn);

    /// \brief PRN for the satellite. Zero if the channel is not allocated.
    int prn = 0;

	double f_carr;                                /*< Carrier frequency */
	double f_code;	/*< Code frequency */
	double carr_phase;
	double code_phase; /*< Code phase */
	gpstime_t g0;	/*!< GPS time at start */
	unsigned long sbf[5][N_DWRD_SBF]; /*!< current subframe */
	unsigned long dwrd[N_DWRD]; /*!< Data words of sub-frame */
	int initial_word;	/*!< initial word */
	int initial_bit;	/*!< initial bit */
	int initial_code;	/*!< initial code */
	AzimuthElevation azel;
	range_t rho0;

    /// \brief Returns true if the channel is enabled.
    bool IsEnabled() const { return prn > 0; }

    int current_data_bit() const { return current_data_bit_; }
    int current_code_chip() const { return current_code_chip_; }

    void UpdateCodeChip() {
        current_code_chip_ = ComputeCodeChip();
    }

    void UpdateDataBit() {
        current_data_bit_ = ComputeDataBit();
    }

    void UpdateCarrierPhase(const double delta_t_s) {
        carr_phase += f_carr * delta_t_s;
        if (carr_phase >= 1.0)
            carr_phase -= 1.0;
        else if (carr_phase<0.0)
            carr_phase += 1.0;
    }

    void UpdateCodePhase(const double delta_t_s)
    {
        code_phase += f_code * delta_t_s;

        if (code_phase >= CA_SEQ_LEN) {
            code_phase -= CA_SEQ_LEN;

            initial_code++;
            if (initial_code>=20) // 20 C/A codes = 1 navigation data bit
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

  private:
    /// \brief C/A code sequence. All values are either -1 or 1.
    std::array<int, CA_SEQ_LEN> code_sequence_;
    
    int current_data_bit_ = 0;
    int current_code_chip_ = 0;
    
    int ComputeDataBit() const {
        return (int)((dwrd[initial_word] >> (29 - initial_bit)) & 0x1UL)*2 - 1;
    }

    int ComputeCodeChip() const {
        return code_sequence_[(int) code_phase];
    }
};

#endif // GPS_CHANNEL_H_

