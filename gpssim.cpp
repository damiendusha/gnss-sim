// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT
//
// Derived from https://github.com/osqzss/gps-sdr-sim (MIT Licence):
// Copyright (c) 2015-2020 Takuji Ebinuma

#include "gpssim.h"
#include "gps_channel.h"
#include "gps_time.h"
#include "gps_ephem.h"
#include "gps_math.h"
#include "gps_subframe.h"
#include "geodesy.h"
#include "ionosphere.h"
#include "noise_generator.h"
#include "rinex2_reader.h"
#include "sin_table.h"
#include "sample_writer.h"
#include "satellite_gain.h"

#include <eigen3/Eigen/Core>

#include <unistd.h>
#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

namespace {
    
static constexpr double kWavelengthGpsL1ca_m = 0.190293672798365;

}   // namespace


/*! \brief Compute range between a satellite and the receiver
 *  \param[out] rho The computed range
 *  \param[in] eph Ephemeris data of the satellite
 *  \param[in] g GPS time at time of receiving the signal
 *  \param[in] e_pos_e_a is the position of the receiver.
 */
void computeRange(range_t *rho, const ephem_t &eph, const ionoutc_t &ionoutc, 
                  gpstime_t g, const Eigen::Vector3d &e_pos_e_a)
{
	double pos[3],vel[3],clk[2];
	double los[3];
	double xrot,yrot;

	double neu[3];
	double tmat[3][3];
	
	// SV position at time of the pseudorange observation.
	satpos(eph, g, pos, vel, clk);

	// Receiver to satellite vector and light-time.
	subVect(los, pos, e_pos_e_a);
	const double tau = normVect(los)/SPEED_OF_LIGHT;

	// Extrapolate the satellite position backwards to the transmission time.
	pos[0] -= vel[0]*tau;
	pos[1] -= vel[1]*tau;
	pos[2] -= vel[2]*tau;

	// Earth rotation correction. The change in velocity can be neglected.
	xrot = pos[0] + pos[1]*OMEGA_EARTH*tau;
	yrot = pos[1] - pos[0]*OMEGA_EARTH*tau;
	pos[0] = xrot;
	pos[1] = yrot;

	// New observer to satellite vector and satellite range.
	subVect(los, pos, e_pos_e_a);
	const double range = normVect(los);
	rho->d = range;

	// Pseudorange.
	rho->range = range - SPEED_OF_LIGHT*clk[0];

	// Relative velocity of SV and receiver.
	const double rate = dotProd(vel, los)/range;

	// Pseudorange rate.
	rho->rate = rate; // - SPEED_OF_LIGHT*clk[1];

	// Time of application.
	rho->g = g;

	// Azimuth and elevation angles.
	const GeodeticPosition llh = xyz2llh(e_pos_e_a);
	ltcmat(llh, tmat);
	ecef2neu(los, tmat, neu);
	rho->azel = neu2azel(neu);

	// Add ionospheric delay
	rho->iono_delay = ionosphericDelay(ionoutc, g, llh, rho->azel);
	rho->range += rho->iono_delay;

	return;
}




void generateNavMsg(gpstime_t g, GpsChannel *chan, int init)
{
	unsigned long previous_word = 0;

    // Align with the full frame length = 30 sec
    gpstime_t g0;
    g0.week = g.week;
    g0.sec = (double)(((unsigned long)(g.sec+0.5))/30UL) * 30.0;
    chan->g0 = g0; // Data bit reference time

    const unsigned long wn = (unsigned long)(g0.week%1024);
    unsigned long tow = ((unsigned long)g0.sec)/6UL;

	if (init == 1) // Initialize subframe 5
	{
		for (int iwrd = 0; iwrd < N_DWRD_SBF; iwrd++)
		{
			unsigned int sbfwrd = chan->sbf[4][iwrd];

			// Add TOW-count message into HOW
			if (iwrd==1)
				sbfwrd |= ((tow&0x1FFFFUL)<<13);

			// Compute checksum
			sbfwrd |= (previous_word<<30) & 0xC0000000UL; // 2 LSBs of the previous transmitted word
			const int nib = ((iwrd==1)||(iwrd==9))?1:0; // Non-information bearing bits for word 2 and 10
			chan->dwrd[iwrd] = computeChecksum(sbfwrd, nib);

			previous_word = chan->dwrd[iwrd];
		}
	}
	else // Save subframe 5
	{
		for (int iwrd = 0; iwrd < N_DWRD_SBF; iwrd++)
		{
			chan->dwrd[iwrd] = chan->dwrd[N_DWRD_SBF*N_SBF+iwrd];

			previous_word = chan->dwrd[iwrd];
		}
	}

	for (int isbf = 0; isbf < N_SBF; isbf++)
	{
		tow++;

		for (int iwrd = 0; iwrd < N_DWRD_SBF; iwrd++)
		{
			unsigned int sbfwrd = chan->sbf[isbf][iwrd];

			// Add transmission week number to Subframe 1
			if ((isbf==0)&&(iwrd==2))
				sbfwrd |= (wn&0x3FFUL)<<20;

			// Add TOW-count message into HOW
			if (iwrd==1)
				sbfwrd |= ((tow&0x1FFFFUL)<<13);

			// Compute checksum
			sbfwrd |= (previous_word<<30) & 0xC0000000UL; // 2 LSBs of the previous transmitted word
			const int nib = ((iwrd==1)||(iwrd==9))?1:0; // Non-information bearing bits for word 2 and 10
			chan->dwrd[(isbf+1)*N_DWRD_SBF+iwrd] = computeChecksum(sbfwrd, nib);

			previous_word = chan->dwrd[(isbf+1)*N_DWRD_SBF+iwrd];
		}
	}
}

bool checkSatVisibility(const ephem_t &eph, gpstime_t g, 
                        const Eigen::Vector3d &e_pos_e_a, 
                        double elevation_mask_deg, AzimuthElevation &out_azel)
{
	double neu[3];
	double pos[3],vel[3],clk[3],los[3];
	double tmat[3][3];

	if (!eph.valid)
        return false;

	const GeodeticPosition llh = xyz2llh(e_pos_e_a);
	ltcmat(llh, tmat);

	satpos(eph, g, pos, vel, clk);
	subVect(los, pos, e_pos_e_a);
	ecef2neu(los, tmat, neu);
	out_azel = neu2azel(neu);

	return out_azel.elevation_deg() > elevation_mask_deg;
}

int allocateChannel(GpsChannel *chan, ephem_t *eph, int* allocatedSat, 
                    ionoutc_t ionoutc, gpstime_t current_simulation_time, 
                    const Eigen::Vector3d &e_pos_e_a, double elevation_mask_deg)
{
	int num_visible_sats = 0;
	AzimuthElevation azel;

	for (int sv = 0; sv < MAX_SAT; sv++)
	{
		if (checkSatVisibility(eph[sv], current_simulation_time, e_pos_e_a, elevation_mask_deg, azel))
		{
			num_visible_sats++;

			if (allocatedSat[sv] == -1) // Visible but not allocated
			{
				// Allocated new satellite
                int i;
				for (i = 0; i < MAX_CHAN; i++)
				{
					if (!chan[i].IsEnabled())
					{
						// Initialize channel
						chan[i] = GpsChannel(sv+1);
						chan[i].azel = azel;

						// Generate subframe
						eph2sbf(eph[sv], ionoutc, chan[i].sbf);

						// Generate navigation message
						generateNavMsg(current_simulation_time, &chan[i], 1);

						// Initialize pseudorange
                        range_t rho;
						computeRange(&rho, eph[sv], ionoutc, 
                                     current_simulation_time, e_pos_e_a);
						chan[i].rho0 = rho;

						// Initialize carrier phase
						const double r_xyz = rho.range;

						computeRange(&rho, eph[sv], ionoutc, 
                                     current_simulation_time, 
                                     Eigen::Vector3d::Zero());
						const double r_ref = rho.range;

						const double phase_ini = 
                                (2.0*r_ref - r_xyz) * (1.0 / kWavelengthGpsL1ca_m);
						chan[i].carrier_phase_cycles = phase_ini - floor(phase_ini);

						// Done.
						break;
					}
				}

				// Set satellite allocation channel
				if (i < MAX_CHAN)
					allocatedSat[sv] = i;
			}
		}
		else if (allocatedSat[sv] >= 0) // Not visible but allocated
		{
			// Clear channel
			chan[allocatedSat[sv]].prn = 0;

			// Clear satellite allocation flag
			allocatedSat[sv] = -1;
		}
	}

	return num_visible_sats;
}

void usage(void)
{
	fprintf(stderr, "Usage: gps-sdr-sim [options]\n"
		"Options:\n"
		"  -e <gps_nav>     RINEX navigation file for GPS ephemerides (required)\n"
		"  -c <location>    ECEF X,Y,Z in meters (static mode) e.g. 3967283.154,1022538.181,4872414.484\n"
		"  -l <location>    Lat,Lon,Hgt (static mode) e.g. 35.681298,139.766247,10.0\n"
		"  -t <date,time>   Scenario start time YYYY/MM/DD,hh:mm:ss\n"
		"  -d <duration>    Duration [sec] (static mode max: %d)\n"
		"  -o <output>      I/Q sampling data file (default: gpssim.bin)\n"
		"  -s <frequency>   Sampling frequency [Hz] (default: 2600000)\n"
		"  -i               Disable ionospheric delay for spacecraft scenario\n"
		"  -v               Show details about simulated channels\n",
		STATIC_MAX_DURATION);

	return;
}

/// \brief Checks whether the sample frequency has an integer number of ns ticks.
bool SampleFrequencyIsValid(int sample_frequency) {
    constexpr int ticks_per_ns = 1000000000;
    int period_ns = ticks_per_ns / sample_frequency;
    return period_ns * sample_frequency == ticks_per_ns;
}

int main(int argc, char *argv[])
{
	clock_t tstart,tend;

	ephem_t eph[EPHEM_ARRAY_SIZE][MAX_SAT];

	gpstime_t simulation_start_gps_time;
    gpstime_t current_simulation_time;
	
    // Default static location; Tokyo
	GeodeticPosition llh = GeodeticPosition::FromDegrees(35.681298, 139.766247, 10.0);

	GpsChannel chan[MAX_CHAN];
	double elevation_mask_deg = 0.0;

    // User position in ECEF coordinates.
    Eigen::Vector3d e_pos_e_a;

    std::string rinex2_nav_file;
    std::string output_sample_filename = {"gpssim.bin"};

	int result;

	datetime_t t0,tmin,tmax;
	gpstime_t gmin,gmax;

	bool verbose = false;

	ionoutc_t ionoutc;
    
    int allocatedSat[MAX_SAT];

    NoiseGenerator noise_generator(1.0, 0.20, 0.5);

	////////////////////////////////////////////////////////////
	// Read options
	////////////////////////////////////////////////////////////

	// Default options
	double raw_samp_freq = 2.6e6;
	simulation_start_gps_time.week = -1; // Invalid start time
	double duration = 300;     // 5 minutes.
	ionoutc.enable = true;

	if (argc<3)
	{
		usage();
		exit(1);
	}

	while ((result=getopt(argc,argv,"e:u:g:c:l:o:s:t:d:iv"))!=-1)
	{
		switch (result)
		{
		case 'e':
            rinex2_nav_file.assign(optarg);
			break;
		case 'c':
        {
			// Static ECEF coordinates input mode
            double xyz[3];
			sscanf(optarg,"%lf,%lf,%lf", &xyz[0], &xyz[1], &xyz[2]);
            e_pos_e_a = {xyz[0], xyz[1], xyz[2]};
			break;
        }
		case 'l':
        {
			// Static geodetic coordinates input mode.
            double raw_lat, raw_lon, raw_height;
			sscanf(optarg,"%lf,%lf,%lf",&raw_lat, &raw_lon, &raw_height);
            llh = GeodeticPosition::FromDegrees(raw_lat, raw_lon, raw_height);
			e_pos_e_a = llh2xyz(llh);
			break;
        }
		case 'o':
			output_sample_filename.assign(optarg);
			break;
		case 's':
			raw_samp_freq = atof(optarg);
			if (raw_samp_freq < 1.0e6)
			{
				fprintf(stderr, "ERROR: Invalid sampling frequency.\n");
				exit(1);
			}
			break;
		case 't':
			sscanf(optarg, "%d/%d/%d,%d:%d:%lf", &t0.y, &t0.m, &t0.d, &t0.hh, &t0.mm, &t0.sec);
			if (t0.y<=1980 || t0.m<1 || t0.m>12 || t0.d<1 || t0.d>31 ||
				t0.hh<0 || t0.hh>23 || t0.mm<0 || t0.mm>59 || t0.sec<0.0 || t0.sec>=60.0)
			{
				fprintf(stderr, "ERROR: Invalid date and time.\n");
				exit(1);
			}
			t0.sec = floor(t0.sec);
			simulation_start_gps_time = date2gps(t0);
			break;
		case 'd':
			duration = atof(optarg);
			break;
		case 'i':
			ionoutc.enable = false; // Disable ionospheric correction
			break;
		case 'v':
			verbose = true;
			break;
		case ':':
		case '?':
			usage();
			exit(1);
		default:
			break;
		}
	}

	if (rinex2_nav_file.empty())
	{
		fprintf(stderr, "ERROR: GPS ephemeris file is not specified.\n");
		exit(1);
	}

	if (duration < 0.0 || duration > STATIC_MAX_DURATION)
	{
		fprintf(stderr, "ERROR: Invalid duration.\n");
		exit(1);
	}

    SampleWriter sample_writer;
    if (sample_writer.OpenFile(output_sample_filename))
    {
        fprintf(stderr, "ERROR: Failed to open output sample file \"%s\".\n", output_sample_filename.c_str());
    }

    const int sample_freq_hz = std::lround(raw_samp_freq);
    if (!SampleFrequencyIsValid(sample_freq_hz)) {
        fprintf(stderr, "ERROR: Sample frequency %f not an integer divisor of 1GHz\n",
                raw_samp_freq);
        exit(1);
    }

    const double sample_period_s = 1.0 / raw_samp_freq;
    const int sample_period_ns = 1000000000 / sample_freq_hz;

    // Currently, samples are processed in batches of 100ms
    const int num_sample_batches = std::lround(duration * 10.0);
    const int samples_per_batch = 100000000 / sample_period_ns;
    const double batch_period_s = 0.1;

	////////////////////////////////////////////////////////////
	// Read ephemeris
	////////////////////////////////////////////////////////////

	const int neph = readRinexNavAll(eph, &ionoutc, rinex2_nav_file);

	if (neph==0)
	{
		fprintf(stderr, "ERROR: No ephemeris available.\n");
		exit(1);
	}

	if (verbose && ionoutc.valid)
	{
		fprintf(stderr, "  %12.3e %12.3e %12.3e %12.3e\n", 
			ionoutc.alpha0, ionoutc.alpha1, ionoutc.alpha2, ionoutc.alpha3);
		fprintf(stderr, "  %12.3e %12.3e %12.3e %12.3e\n", 
			ionoutc.beta0, ionoutc.beta1, ionoutc.beta2, ionoutc.beta3);
		fprintf(stderr, "   %19.11e %19.11e  %9d %9d\n",
			ionoutc.A0, ionoutc.A1, ionoutc.tot, ionoutc.wnt);
		fprintf(stderr, "%6d\n", ionoutc.dtls);
	}

	for (int sv = 0; sv < MAX_SAT; sv++) 
	{
		if (eph[0][sv].valid)
		{
			gmin = eph[0][sv].toc;
			tmin = eph[0][sv].t;
			break;
		}
	}

	gmax.sec = 0;
	gmax.week = 0;
	tmax.sec = 0;
	tmax.mm = 0;
	tmax.hh = 0;
	tmax.d = 0;
	tmax.m = 0;
	tmax.y = 0;
	for (int sv = 0; sv < MAX_SAT; sv++)
	{
		if (eph[neph-1][sv].valid)
		{
			gmax = eph[neph-1][sv].toc;
			tmax = eph[neph-1][sv].t;
			break;
		}
	}

	// Scenario start time has been set.
	if (simulation_start_gps_time.week >= 0)
	{
		if (subGpsTime(simulation_start_gps_time, gmin) < 0.0 || 
            subGpsTime(gmax, simulation_start_gps_time) < 0.0)
		{
			fprintf(stderr, "ERROR: Invalid start time.\n");
			fprintf(stderr, "tmin = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n", 
				tmin.y, tmin.m, tmin.d, tmin.hh, tmin.mm, tmin.sec,
				gmin.week, gmin.sec);
			fprintf(stderr, "tmax = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n", 
				tmax.y, tmax.m, tmax.d, tmax.hh, tmax.mm, tmax.sec,
				gmax.week, gmax.sec);
			exit(1);
		}
	}
	else
	{
		simulation_start_gps_time = gmin;
		t0 = tmin;
	}

	fprintf(stderr, "Start time = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n", 
		t0.y, t0.m, t0.d, t0.hh, t0.mm, t0.sec, simulation_start_gps_time.week, 
         simulation_start_gps_time.sec);
	fprintf(stderr, "Duration = %.1f [sec]\n", duration);

	// Select the current set of ephemerides
	int ieph = -1;

	for (int i = 0; i < neph; i++)
	{
		for (int sv = 0; sv < MAX_SAT; sv++)
		{
			if (eph[i][sv].valid)
			{
				const double dt = subGpsTime(simulation_start_gps_time, eph[i][sv].toc);
				if (dt>=-SECONDS_IN_HOUR && dt<SECONDS_IN_HOUR)
				{
					ieph = i;
					break;
				}
			}
		}

		if (ieph>=0) // ieph has been set
			break;
	}

	if (ieph == -1)
	{
		fprintf(stderr, "ERROR: No current set of ephemerides has been found.\n");
		exit(1);
	}

	////////////////////////////////////////////////////////////
	// Initialize channels
	////////////////////////////////////////////////////////////

	// Clear satellite allocation flag
	for (int sv = 0; sv < MAX_SAT; sv++) {
		allocatedSat[sv] = -1;
    }

	// Initial reception time
	current_simulation_time = incGpsTime(simulation_start_gps_time, 0.0);

	// Allocate visible satellites
	allocateChannel(chan, eph[ieph], allocatedSat, ionoutc, current_simulation_time, 
                    e_pos_e_a, elevation_mask_deg);

	for(int i = 0; i < MAX_CHAN; i++)
	{
		if (chan[i].IsEnabled()) {
			fprintf(stderr, "%02d %6.1f %5.1f %11.1f %5.1f\n", chan[i].prn, 
				chan[i].azel.azimuth_deg(), chan[i].azel.elevation_deg(), chan[i].
				rho0.d, chan[i].rho0.iono_delay);
        }
	}

	////////////////////////////////////////////////////////////
	// Generate baseband signals
	////////////////////////////////////////////////////////////

	tstart = clock();

	// Update receiver time
	current_simulation_time = incGpsTime(current_simulation_time, batch_period_s);
    
    PeriodicFunctionTable<double, 1024> sin_table(
        [](double t) -> double {return std::sin(t);}, 1.0);
    PeriodicFunctionTable<double, 1024> cos_table(
        [](double t) -> double {return std::cos(t);}, 1.0);
    
    // Default all channels to a gain of 1.
    ConstellationGain constellation_gain;
    for (int prn = 0; prn <= 32; ++prn) {
        constellation_gain.SetSatelliteToConstantGain(prn, 1.0);
    }
    
    // Successively decrease the gain as the PRN increases.
    constellation_gain.SetSatelliteToConstantGain( 1, 1.0);
    constellation_gain.SetSatelliteToConstantGain( 2, 0.9);
    constellation_gain.SetSatelliteToConstantGain( 3, 0.7);
    constellation_gain.SetSatelliteToConstantGain( 6, 0.6);
    constellation_gain.SetSatelliteToConstantGain( 9, 0.5);
    constellation_gain.SetSatelliteToConstantGain(10, 0.4);
    constellation_gain.SetSatelliteToConstantGain(11, 0.3);
    constellation_gain.SetSatelliteToConstantGain(12, 0.2);
    constellation_gain.SetSatelliteToConstantGain(17, 0.17);
    constellation_gain.SetSatelliteToConstantGain(20, 0.15);
    constellation_gain.SetSatelliteToConstantGain(23, 0.13);
    constellation_gain.SetSatelliteToConstantGain(28, 0.10);
    constellation_gain.SetSatelliteToConstantGain(32, 0.07);

	for (int batch = 1; batch < num_sample_batches; batch++)
	{
        // Per-channel gain, maximum value of 1.
        std::array<double, MAX_CHAN> gain;

		for (int i = 0; i < MAX_CHAN; i++)
		{
			if (chan[i].IsEnabled())
			{
				// Refresh code phase and data bit counters
				range_t rho;
				const int sv_index = chan[i].prn-1;

				// Current pseudorange
				computeRange(&rho, eph[ieph][sv_index], ionoutc, 
                             current_simulation_time, e_pos_e_a);

				chan[i].azel = rho.azel;

				// Update code phase and data bit counters
				chan[i].ComputeCodePhase(rho, batch_period_s);

				// Signal gain.
				gain[i] = constellation_gain.ComputeGain(
                    chan[i].prn, rho.d, rho.azel);
			}
		}

		// There is no intrinsic need to loop here, but we do so to avoid computing
		// the pseudorange too frequently.
		for (int sample = 0; sample < samples_per_batch; sample++)
		{
			double i_acc = 0;
			double q_acc = 0;

			for (int i = 0; i < MAX_CHAN; i++)
			{
				if (chan[i].IsEnabled())
				{
                    const double coeff = chan[i].current_data_bit() * 
                            chan[i].current_code_chip() * gain[i];
                    const double ip = coeff * 
                            cos_table.LookupValue(chan[i].carrier_phase_cycles);
                    const double qp = coeff * 
                            sin_table.LookupValue(chan[i].carrier_phase_cycles);

					// Accumulate for all visible satellites.
					i_acc += ip;
					q_acc += qp;

					// Update code phase, including the code chip and data bit.
					chan[i].UpdateCodePhase(sample_period_s);

					// Update carrier phase
					chan[i].UpdateCarrierPhase(sample_period_s);
				}
			}

            // At this point, the min and max values for the I- and Q- samples
            // are -16 and +16 based on the number of satellites processed.
            //
            // We scale them to +/- 0.5 for the sample writer.
            constexpr double kSampleScale = 1.0 / (2.0 * 16.0);
            sample_writer.WriteSample(
                    noise_generator.ScaleAndAddNoise(i_acc * kSampleScale),
                    noise_generator.ScaleAndAddNoise(q_acc * kSampleScale));
        }

        //
        // Update navigation message and channel allocation on 30s boundaries.
        //

        const int icurrent_simulation_time = 
            (int)(current_simulation_time.sec*10.0+0.5);

        if (icurrent_simulation_time % 300 == 0)
		{
			// Update navigation message
			for (int i = 0; i < MAX_CHAN; i++)
			{
				if (chan[i].prn>0)
					generateNavMsg(current_simulation_time, &chan[i], 0);
			}

			// Refresh ephemeris and subframes
			// Quick and dirty fix. Need more elegant way.
			for (int sv = 0; sv < MAX_SAT; sv++)
			{
				if (eph[ieph+1][sv].valid)
				{
					const double dt = subGpsTime(eph[ieph+1][sv].toc, 
                                                 current_simulation_time);
					if (dt < SECONDS_IN_HOUR)
					{
						ieph++;
						for (int i = 0; i < MAX_CHAN; i++)
						{
							// Generate new subframes if allocated
							if (chan[i].IsEnabled()) { 
								eph2sbf(eph[ieph][chan[i].prn-1], ionoutc, chan[i].sbf);
                            }
						}
					}
						
					break;
				}
			}

			// Update channel allocation
			allocateChannel(chan, eph[ieph], allocatedSat, ionoutc, 
                            current_simulation_time, e_pos_e_a, elevation_mask_deg);

			// Show details about simulated channels
			if (verbose)
			{
				fprintf(stderr, "\n");
				for (int i = 0; i < MAX_CHAN; i++)
				{
					if (chan[i].prn>0)
						fprintf(stderr, "%02d %6.1f %5.1f %11.1f %5.1f\n", chan[i].prn,
							chan[i].azel.azimuth_deg(), chan[i].azel.elevation_deg(),
                            chan[i].rho0.d, chan[i].rho0.iono_delay);
				}
			}
		}

		// Update receiver time
		current_simulation_time = incGpsTime(current_simulation_time, batch_period_s);

		// Update time counter
		fprintf(stderr, "\rTime into run = %4.1f", subGpsTime(current_simulation_time, simulation_start_gps_time));
		fflush(stdout);
	}

	tend = clock();

	fprintf(stderr, "\nDone!\n");

	// Process time
	fprintf(stderr, "Process time = %.1f [sec]\n", (double)(tend-tstart)/CLOCKS_PER_SEC);

	return(0);
}
