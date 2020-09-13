// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT
//
// Derived from https://github.com/osqzss/gps-sdr-sim (MIT Licence):
// Copyright (c) 2015-2020 Takuji Ebinuma

#include "gps_channel.h"
#include "gps_subframe.h"

#include <cassert>


namespace {
    
static constexpr double kWavelengthGpsL1ca_m = 0.190293672798365;

/*! \brief generate the C/A code sequence for a given Satellite Vehicle PRN
 *  \param[in] prn PRN nuber of the Satellite Vehicle
 *  \param[out] ca Caller-allocated integer array of 1023 bytes
 */
std::array<int, CA_SEQ_LEN> GeneratePrnSequence(int prn)
{
    assert(prn >= 1 && prn <= 32);
    std::array<int, CA_SEQ_LEN> ca;

	constexpr std::array<int, 32> delay = {
		  5,   6,   7,   8,  17,  18, 139, 140, 141, 251,
		252, 254, 255, 256, 257, 258, 469, 470, 471, 472,
		473, 474, 509, 512, 513, 514, 515, 516, 859, 860,
		861, 862};
	
	int g1[CA_SEQ_LEN], g2[CA_SEQ_LEN];
	int r1[N_DWRD_SBF], r2[N_DWRD_SBF];

	for (int i = 0; i < N_DWRD_SBF; i++) {
		r1[i] = r2[i] = -1;
    }

	for (int i = 0; i < CA_SEQ_LEN; i++)
	{
		g1[i] = r1[9];
		g2[i] = r2[9];
		int c1 = r1[2]*r1[9];
		int c2 = r2[1]*r2[2]*r2[5]*r2[7]*r2[8]*r2[9];

		for (int j = 9; j > 0; j--) 
		{
			r1[j] = r1[j-1];
			r2[j] = r2[j-1];
		}
		r1[0] = c1;
		r2[0] = c2;
	}

	for (int i = 0, j = CA_SEQ_LEN-delay[prn-1]; i < CA_SEQ_LEN; i++,j++) {
		ca[i] = (1-g1[i]*g2[j%CA_SEQ_LEN])/2 == 0 ? -1 : 1;
    }

	return ca;
}

}   // namespace

GpsChannel::GpsChannel(int in_prn)
    : prn(in_prn)
    , code_sequence_(GeneratePrnSequence(in_prn))
{
}

/*! \brief Compute the code phase for a given channel (satellite)
 *  \param chan Channel on which we operate (is updated)
 *  \param[in] rho1 Current range, after \a dt has expired
 *  \param[in dt delta-t (time difference) in seconds
 */
void GpsChannel::ComputeCodePhase(const range_t &rho1, const double dt)
{
	// Pseudorange rate.
    // TODO: The range rate is computed in range_t, but not applied here
    // Why is that?
	const double rhorate = (rho1.range - rho0.range) / dt;

	// Carrier and code frequency.
	f_carr = -rhorate * (1.0 / kWavelengthGpsL1ca_m);
	f_code = CODE_FREQ + f_carr*CARR_TO_CODE;

	// Initial code phase and data bit counters.
	const double ms = ((subGpsTime(rho0.g, dataframe_reference_time)+6.0) - 
            rho0.range/SPEED_OF_LIGHT)*1000.0;

	int ims = static_cast<int>(std::floor(ms));
	code_phase = (ms-(double)ims)*CA_SEQ_LEN; // in chip

	initial_word = ims/600; // 1 word = 30 bits = 600 ms
	ims -= initial_word*600;

	initial_bit = ims/20; // 1 bit = 20 code = 20 ms
	ims -= initial_bit*20;

	initial_code = ims; // 1 code = 1 ms

	UpdateCodeChip();
	UpdateDataBit();

	// Save current pseudorange
	rho0 = rho1;
}

void GpsChannel::GenerateNavMsg(gpstime_t g, int init)
{
	unsigned long previous_word = 0;

    // Align with the full frame length = 30 sec
    gpstime_t g0;
    g0.week = g.week;
    g0.sec = (double)(((unsigned long)(g.sec+0.5))/30UL) * 30.0;

    dataframe_reference_time = g0; // Data bit reference time

    const unsigned long wn = (unsigned long)(g0.week%1024);
    unsigned long tow = ((unsigned long)g0.sec)/6UL;

	if (init == 1) // Initialize subframe 5
	{
		for (int iwrd = 0; iwrd < N_DWRD_SBF; iwrd++)
		{
			unsigned int sbfwrd = sbf[4][iwrd];

			// Add TOW-count message into HOW
			if (iwrd==1)
				sbfwrd |= ((tow&0x1FFFFUL)<<13);

			// Compute checksum
			sbfwrd |= (previous_word<<30) & 0xC0000000UL; // 2 LSBs of the previous transmitted word
			const int nib = ((iwrd==1)||(iwrd==9))?1:0; // Non-information bearing bits for word 2 and 10
			dwrd[iwrd] = computeChecksum(sbfwrd, nib);

			previous_word = dwrd[iwrd];
		}
	}
	else // Save subframe 5
	{
		for (int iwrd = 0; iwrd < N_DWRD_SBF; iwrd++)
		{
			dwrd[iwrd] = dwrd[N_DWRD_SBF*N_SBF+iwrd];

			previous_word = dwrd[iwrd];
		}
	}

	for (int isbf = 0; isbf < N_SBF; isbf++)
	{
		tow++;

		for (int iwrd = 0; iwrd < N_DWRD_SBF; iwrd++)
		{
			unsigned int sbfwrd = sbf[isbf][iwrd];

			// Add transmission week number to Subframe 1
			if ((isbf==0)&&(iwrd==2))
				sbfwrd |= (wn&0x3FFUL)<<20;

			// Add TOW-count message into HOW
			if (iwrd==1)
				sbfwrd |= ((tow&0x1FFFFUL)<<13);

			// Compute checksum
			sbfwrd |= (previous_word<<30) & 0xC0000000UL; // 2 LSBs of the previous transmitted word
			const int nib = ((iwrd==1)||(iwrd==9))?1:0; // Non-information bearing bits for word 2 and 10
			dwrd[(isbf+1)*N_DWRD_SBF+iwrd] = computeChecksum(sbfwrd, nib);

			previous_word = dwrd[(isbf+1)*N_DWRD_SBF+iwrd];
		}
	}
}

