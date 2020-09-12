// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT
//
// Derived from https://github.com/osqzss/gps-sdr-sim (MIT Licence):
// Copyright (c) 2015-2020 Takuji Ebinuma

#include "gpssim.h"
#include "rinex2_reader.h"
#include "gps_time.h"
#include "gps_ephem.h"
#include "gps_math.h"
#include "gps_subframe.h"
#include "geodesy.h"
#include "noise_generator.h"
#include "sin_table.h"
#include "sample_writer.h"
#include "satellite_gain.h"

#include <unistd.h>
#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>


int allocatedSat[MAX_SAT];


/* !\brief generate the C/A code sequence for a given Satellite Vehicle PRN
 *  \param[in] prn PRN nuber of the Satellite Vehicle
 *  \param[out] ca Caller-allocated integer array of 1023 bytes
 */
void codegen(int *ca, int prn)
{
	int delay[] = {
		  5,   6,   7,   8,  17,  18, 139, 140, 141, 251,
		252, 254, 255, 256, 257, 258, 469, 470, 471, 472,
		473, 474, 509, 512, 513, 514, 515, 516, 859, 860,
		861, 862};
	
	int g1[CA_SEQ_LEN], g2[CA_SEQ_LEN];
	int r1[N_DWRD_SBF], r2[N_DWRD_SBF];

	if (prn < 1 || prn > 32) {
		return;
    }

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
		ca[i] = (1-g1[i]*g2[j%CA_SEQ_LEN])/2;
    }

	return;
}


double ionosphericDelay(const ionoutc_t *ionoutc, gpstime_t g, double *llh, double *azel)
{
	if (ionoutc->enable==false)
		return (0.0); // No ionospheric delay

	const double E = azel[1]/PI;
	const double phi_u = llh[0]/PI;
	const double lam_u = llh[1]/PI;

	// Obliquity factor
	const double F = 1.0 + 16.0*pow((0.53 - E),3.0);

    double iono_delay = 0.0;
	if (ionoutc->vflg==false) {
		iono_delay = F*5.0e-9*SPEED_OF_LIGHT;
    }
	else
	{
		// Earth's central angle between the user position and the earth projection of
		// ionospheric intersection point (semi-circles)
		const double psi = 0.0137/(E + 0.11) - 0.022;
		
		// Geodetic latitude of the earth projection of the ionospheric intersection point
		// (semi-circles)
		double phi_i = phi_u + psi*cos(azel[0]);
		if(phi_i>0.416)
			phi_i = 0.416;
		else if(phi_i<-0.416)
			phi_i = -0.416;

		// Geodetic longitude of the earth projection of the ionospheric intersection point
		// (semi-circles)
		const double lam_i = lam_u + psi*sin(azel[0])/cos(phi_i*PI);

		// Geomagnetic latitude of the earth projection of the ionospheric intersection
		// point (mean ionospheric height assumed 350 km) (semi-circles)
		const double phi_m = phi_i + 0.064*cos((lam_i - 1.617)*PI);
		const double phi_m2 = phi_m*phi_m;
		const double phi_m3 = phi_m2*phi_m;

		double AMP = ionoutc->alpha0 + ionoutc->alpha1*phi_m
			+ ionoutc->alpha2*phi_m2 + ionoutc->alpha3*phi_m3;
		if (AMP<0.0)
			AMP = 0.0;

		double PER = ionoutc->beta0 + ionoutc->beta1*phi_m
			+ ionoutc->beta2*phi_m2 + ionoutc->beta3*phi_m3;
		if (PER<72000.0)
			PER = 72000.0;

		// Local time (sec)
		double t = SECONDS_IN_DAY/2.0*lam_i + g.sec;
		while(t>=SECONDS_IN_DAY)
			t -= SECONDS_IN_DAY;
		while(t<0)
			t += SECONDS_IN_DAY;

		// Phase (radians)
		const double X = 2.0*PI*(t - 50400.0)/PER;

		if(fabs(X)<1.57)
		{
			const double X2 = X*X;
			const double X4 = X2*X2;
			iono_delay = F*(5.0e-9 + AMP*(1.0 - X2/2.0 + X4/24.0))*SPEED_OF_LIGHT;
		}
		else
			iono_delay = F*5.0e-9*SPEED_OF_LIGHT;
	}

	return iono_delay;
}

/*! \brief Compute range between a satellite and the receiver
 *  \param[out] rho The computed range
 *  \param[in] eph Ephemeris data of the satellite
 *  \param[in] g GPS time at time of receiving the signal
 *  \param[in] xyz position of the receiver
 */
void computeRange(range_t *rho, const ephem_t &eph, ionoutc_t *ionoutc, gpstime_t g, double xyz[])
{
	double pos[3],vel[3],clk[2];
	double los[3];
	double range,rate;
	double xrot,yrot;

	double llh[3],neu[3];
	double tmat[3][3];
	
	// SV position at time of the pseudorange observation.
	satpos(eph, g, pos, vel, clk);

	// Receiver to satellite vector and light-time.
	subVect(los, pos, xyz);
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
	subVect(los, pos, xyz);
	range = normVect(los);
	rho->d = range;

	// Pseudorange.
	rho->range = range - SPEED_OF_LIGHT*clk[0];

	// Relative velocity of SV and receiver.
	rate = dotProd(vel, los)/range;

	// Pseudorange rate.
	rho->rate = rate; // - SPEED_OF_LIGHT*clk[1];

	// Time of application.
	rho->g = g;

	// Azimuth and elevation angles.
	xyz2llh(xyz, llh);
	ltcmat(llh, tmat);
	ecef2neu(los, tmat, neu);
	neu2azel(rho->azel, neu);

	// Add ionospheric delay
	rho->iono_delay = ionosphericDelay(ionoutc, g, llh, rho->azel);
	rho->range += rho->iono_delay;

	return;
}

/*! \brief Compute the code phase for a given channel (satellite)
 *  \param chan Channel on which we operate (is updated)
 *  \param[in] rho1 Current range, after \a dt has expired
 *  \param[in dt delta-t (time difference) in seconds
 */
void computeCodePhase(channel_t *chan, range_t rho1, double dt)
{
	// Pseudorange rate.
	const double rhorate = (rho1.range - chan->rho0.range)/dt;

	// Carrier and code frequency.
	chan->f_carr = -rhorate/LAMBDA_L1;
	chan->f_code = CODE_FREQ + chan->f_carr*CARR_TO_CODE;

	// Initial code phase and data bit counters.
	const double ms = ((subGpsTime(chan->rho0.g,chan->g0)+6.0) - chan->rho0.range/SPEED_OF_LIGHT)*1000.0;

	int ims = static_cast<int>(ms);
	chan->code_phase = (ms-(double)ims)*CA_SEQ_LEN; // in chip

	chan->iword = ims/600; // 1 word = 30 bits = 600 ms
	ims -= chan->iword*600;
			
	chan->ibit = ims/20; // 1 bit = 20 code = 20 ms
	ims -= chan->ibit*20;

	chan->icode = ims; // 1 code = 1 ms

	chan->codeCA = chan->ca[(int)chan->code_phase]*2-1;
	chan->dataBit = (int)((chan->dwrd[chan->iword]>>(29-chan->ibit)) & 0x1UL)*2-1;

	// Save current pseudorange
	chan->rho0 = rho1;

	return;
}


int generateNavMsg(gpstime_t g, channel_t *chan, int init)
{
	int iwrd,isbf;
	gpstime_t g0;
	unsigned long wn,tow;
	unsigned sbfwrd;
	unsigned long prevwrd;
	int nib;

	g0.week = g.week;
	g0.sec = (double)(((unsigned long)(g.sec+0.5))/30UL) * 30.0; // Align with the full frame length = 30 sec
	chan->g0 = g0; // Data bit reference time

	wn = (unsigned long)(g0.week%1024);
	tow = ((unsigned long)g0.sec)/6UL;

	if (init==1) // Initialize subframe 5
	{
		prevwrd = 0UL;

		for (iwrd=0; iwrd<N_DWRD_SBF; iwrd++)
		{
			sbfwrd = chan->sbf[4][iwrd];

			// Add TOW-count message into HOW
			if (iwrd==1)
				sbfwrd |= ((tow&0x1FFFFUL)<<13);

			// Compute checksum
			sbfwrd |= (prevwrd<<30) & 0xC0000000UL; // 2 LSBs of the previous transmitted word
			nib = ((iwrd==1)||(iwrd==9))?1:0; // Non-information bearing bits for word 2 and 10
			chan->dwrd[iwrd] = computeChecksum(sbfwrd, nib);

			prevwrd = chan->dwrd[iwrd];
		}
	}
	else // Save subframe 5
	{
		for (iwrd=0; iwrd<N_DWRD_SBF; iwrd++)
		{
			chan->dwrd[iwrd] = chan->dwrd[N_DWRD_SBF*N_SBF+iwrd];

			prevwrd = chan->dwrd[iwrd];
		}
		/*
		// Sanity check
		if (((chan->dwrd[1])&(0x1FFFFUL<<13)) != ((tow&0x1FFFFUL)<<13))
		{
			fprintf(stderr, "\nWARNING: Invalid TOW in subframe 5.\n");
			return(0);
		}
		*/
	}

	for (isbf=0; isbf<N_SBF; isbf++)
	{
		tow++;

		for (iwrd=0; iwrd<N_DWRD_SBF; iwrd++)
		{
			sbfwrd = chan->sbf[isbf][iwrd];

			// Add transmission week number to Subframe 1
			if ((isbf==0)&&(iwrd==2))
				sbfwrd |= (wn&0x3FFUL)<<20;

			// Add TOW-count message into HOW
			if (iwrd==1)
				sbfwrd |= ((tow&0x1FFFFUL)<<13);

			// Compute checksum
			sbfwrd |= (prevwrd<<30) & 0xC0000000UL; // 2 LSBs of the previous transmitted word
			nib = ((iwrd==1)||(iwrd==9))?1:0; // Non-information bearing bits for word 2 and 10
			chan->dwrd[(isbf+1)*N_DWRD_SBF+iwrd] = computeChecksum(sbfwrd, nib);

			prevwrd = chan->dwrd[(isbf+1)*N_DWRD_SBF+iwrd];
		}
	}

	return(1);
}

bool checkSatVisibility(ephem_t eph, gpstime_t g, double *xyz, double elevation_mask_deg, double *azel)
{
	double llh[3],neu[3];
	double pos[3],vel[3],clk[3],los[3];
	double tmat[3][3];

	if (eph.vflg != 1) return false;

	xyz2llh(xyz,llh);
	ltcmat(llh, tmat);

	satpos(eph, g, pos, vel, clk);
	subVect(los, pos, xyz);
	ecef2neu(los, tmat, neu);
	neu2azel(azel, neu);

	return azel[1] * R2D > elevation_mask_deg;
}

int allocateChannel(channel_t *chan, ephem_t *eph, ionoutc_t ionoutc, 
                    gpstime_t grx, double *xyz, double elevation_mask_deg)
{
	int num_visible_sats = 0;
	double azel[2];

	double ref[3]={0.0};

	for (int sv = 0; sv < MAX_SAT; sv++)
	{
		if (checkSatVisibility(eph[sv], grx, xyz, elevation_mask_deg, azel))
		{
			num_visible_sats++;

			if (allocatedSat[sv]==-1) // Visible but not allocated
			{
				// Allocated new satellite
                int i;
				for (i = 0; i < MAX_CHAN; i++)
				{
					if (chan[i].prn==0)
					{
						// Initialize channel
						chan[i].prn = sv+1;
						chan[i].azel[0] = azel[0];
						chan[i].azel[1] = azel[1];

						// C/A code generation
						codegen(chan[i].ca, chan[i].prn);

						// Generate subframe
						eph2sbf(eph[sv], ionoutc, chan[i].sbf);

						// Generate navigation message
						generateNavMsg(grx, &chan[i], 1);

						// Initialize pseudorange
                        range_t rho;
						computeRange(&rho, eph[sv], &ionoutc, grx, xyz);
						chan[i].rho0 = rho;

						// Initialize carrier phase
						const double r_xyz = rho.range;

						computeRange(&rho, eph[sv], &ionoutc, grx, ref);
						const double r_ref = rho.range;

						const double phase_ini = (2.0*r_ref - r_xyz)/LAMBDA_L1;
						chan[i].carr_phase = phase_ini - floor(phase_ini);

						// Done.
						break;
					}
				}

				// Set satellite allocation channel
				if (i < MAX_CHAN)
					allocatedSat[sv] = i;
			}
		}
		else if (allocatedSat[sv]>=0) // Not visible but allocated
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
		"  -T <date,time>   Overwrite TOC and TOE to scenario start time\n"
		"  -d <duration>    Duration [sec] (dynamic mode max: %.0f, static mode max: %d)\n"
		"  -o <output>      I/Q sampling data file (default: gpssim.bin)\n"
		"  -s <frequency>   Sampling frequency [Hz] (default: 2600000)\n"
		"  -i               Disable ionospheric delay for spacecraft scenario\n"
		"  -v               Show details about simulated channels\n",
		((double)USER_MOTION_SIZE) / 10.0, STATIC_MAX_DURATION);

	return;
}

int main(int argc, char *argv[])
{
	clock_t tstart,tend;

	int neph,ieph;
	ephem_t eph[EPHEM_ARRAY_SIZE][MAX_SAT];
	gpstime_t g0;
	
	double llh[3];

	channel_t chan[MAX_CHAN];
	double elevation_mask_deg = 0.0;

	gpstime_t grx;

	char umfile[MAX_CHAR];
	double xyz[USER_MOTION_SIZE][3];

	char navfile[MAX_CHAR];

    std::string output_sample_filename = {"gpssim.bin"};

	int result;

	datetime_t t0,tmin,tmax;
	gpstime_t gmin,gmax;

	bool verbose = false;

    // Overwrite the TOC and TOE in the RINEX file.
	bool timeoverwrite = false;

	ionoutc_t ionoutc;
    
    NoiseGenerator noise_generator(1.0, 0.20, 0.5);

	////////////////////////////////////////////////////////////
	// Read options
	////////////////////////////////////////////////////////////

	// Default options
	navfile[0] = '\0';
	umfile[0] = '\0';
	double samp_freq = 2.6e6;
	g0.week = -1; // Invalid start time
	double duration = ((double) USER_MOTION_SIZE)/10.0; // Default duration
	ionoutc.enable = true;

	if (argc<3)
	{
		usage();
		exit(1);
	}

	while ((result=getopt(argc,argv,"e:u:g:c:l:o:s:T:t:d:iv"))!=-1)
	{
		switch (result)
		{
		case 'e':
			strcpy(navfile, optarg);
			break;
		case 'c':
			// Static ECEF coordinates input mode
			sscanf(optarg,"%lf,%lf,%lf",&xyz[0][0],&xyz[0][1],&xyz[0][2]);
			break;
		case 'l':
			// Static geodetic coordinates input mode
			sscanf(optarg,"%lf,%lf,%lf",&llh[0],&llh[1],&llh[2]);
			llh[0] = llh[0] / R2D; // convert to RAD
			llh[1] = llh[1] / R2D; // convert to RAD
			llh2xyz(llh,xyz[0]); // Convert llh to xyz
			break;
		case 'o':
			output_sample_filename.assign(optarg);
			break;
		case 's':
			samp_freq = atof(optarg);
			if (samp_freq < 1.0e6)
			{
				fprintf(stderr, "ERROR: Invalid sampling frequency.\n");
				exit(1);
			}
			break;
		case 'T':
			timeoverwrite = true;
			if (strncmp(optarg, "now", 3)==0)
			{
				time_t timer;
				struct tm *gmt;
				
				time(&timer);
				gmt = gmtime(&timer);

				t0.y = gmt->tm_year+1900;
				t0.m = gmt->tm_mon+1;
				t0.d = gmt->tm_mday;
				t0.hh = gmt->tm_hour;
				t0.mm = gmt->tm_min;
				t0.sec = (double)gmt->tm_sec;

				date2gps(&t0, &g0);

				break;
			}
		case 't':
			sscanf(optarg, "%d/%d/%d,%d:%d:%lf", &t0.y, &t0.m, &t0.d, &t0.hh, &t0.mm, &t0.sec);
			if (t0.y<=1980 || t0.m<1 || t0.m>12 || t0.d<1 || t0.d>31 ||
				t0.hh<0 || t0.hh>23 || t0.mm<0 || t0.mm>59 || t0.sec<0.0 || t0.sec>=60.0)
			{
				fprintf(stderr, "ERROR: Invalid date and time.\n");
				exit(1);
			}
			t0.sec = floor(t0.sec);
			date2gps(&t0, &g0);
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

	if (navfile[0]==0)
	{
		fprintf(stderr, "ERROR: GPS ephemeris file is not specified.\n");
		exit(1);
	}

	// Default static location; Tokyo
	llh[0] = 35.681298 / R2D;
	llh[1] = 139.766247 / R2D;
	llh[2] = 10.0;

	if (duration < 0.0 || duration > STATIC_MAX_DURATION)
	{
		fprintf(stderr, "ERROR: Invalid duration.\n");
		exit(1);
	}
	const int iduration = (int)(duration*10.0 + 0.5);

    SampleWriter sample_writer;
    if (sample_writer.OpenFile(output_sample_filename))
    {
        fprintf(stderr, "ERROR: Failed to open output sample file \"%s\".\n", output_sample_filename.c_str());
    }

	// Round the sample frequency to the nearest 10Hz.
	samp_freq = std::round(samp_freq / 10.0) * 10.0;
	const double delt = 1.0 / samp_freq;

	////////////////////////////////////////////////////////////
	// Receiver position
	////////////////////////////////////////////////////////////

    // Always use a static location.
    const int numd = iduration;

	////////////////////////////////////////////////////////////
	// Read ephemeris
	////////////////////////////////////////////////////////////

	neph = readRinexNavAll(eph, &ionoutc, navfile);

	if (neph==0)
	{
		fprintf(stderr, "ERROR: No ephemeris available.\n");
		exit(1);
	}

	if (verbose && (ionoutc.vflg==true))
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
		if (eph[0][sv].vflg==1)
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
		if (eph[neph-1][sv].vflg == 1)
		{
			gmax = eph[neph-1][sv].toc;
			tmax = eph[neph-1][sv].t;
			break;
		}
	}

	if (g0.week>=0) // Scenario start time has been set.
	{
		if (timeoverwrite==true)
		{
			gpstime_t gtmp;
			datetime_t ttmp;
			double dsec;

			gtmp.week = g0.week;
			gtmp.sec = (double)(((int)(g0.sec))/7200)*7200.0;

			dsec = subGpsTime(gtmp,gmin);

			// Overwrite the UTC reference week number
			ionoutc.wnt = gtmp.week;
			ionoutc.tot = (int)gtmp.sec;

			// Overwrite the TOC and TOE to the scenario start time
			for (int sv = 0; sv < MAX_SAT; sv++)
			{
				for (int i = 0; i < neph; i++)
				{
					if (eph[i][sv].vflg == 1)
					{
						gtmp = incGpsTime(eph[i][sv].toc, dsec);
						gps2date(&gtmp,&ttmp);
						eph[i][sv].toc = gtmp;
						eph[i][sv].t = ttmp;

						gtmp = incGpsTime(eph[i][sv].toe, dsec);
						eph[i][sv].toe = gtmp;
					}
				}
			}
		}
		else
		{
			if (subGpsTime(g0, gmin)<0.0 || subGpsTime(gmax, g0)<0.0)
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
	}
	else
	{
		g0 = gmin;
		t0 = tmin;
	}

	fprintf(stderr, "Start time = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n", 
		t0.y, t0.m, t0.d, t0.hh, t0.mm, t0.sec, g0.week, g0.sec);
	fprintf(stderr, "Duration = %.1f [sec]\n", ((double)numd)/10.0);

	// Select the current set of ephemerides
	ieph = -1;

	for (int i = 0; i < neph; i++)
	{
		for (int sv = 0; sv < MAX_SAT; sv++)
		{
			if (eph[i][sv].vflg == 1)
			{
				const double dt = subGpsTime(g0, eph[i][sv].toc);
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

	// Clear all channels
	for (int i = 0; i < MAX_CHAN; i++) {
		chan[i].prn = 0;
    }

	// Clear satellite allocation flag
	for (int sv = 0; sv < MAX_SAT; sv++) {
		allocatedSat[sv] = -1;
    }

	// Initial reception time
	grx = incGpsTime(g0, 0.0);

	// Allocate visible satellites
	allocateChannel(chan, eph[ieph], ionoutc, grx, xyz[0], elevation_mask_deg);

	for(int i = 0; i < MAX_CHAN; i++)
	{
		if (chan[i].prn>0)
			fprintf(stderr, "%02d %6.1f %5.1f %11.1f %5.1f\n", chan[i].prn, 
				chan[i].azel[0]*R2D, chan[i].azel[1]*R2D, chan[i].rho0.d, chan[i].rho0.iono_delay);
	}

	////////////////////////////////////////////////////////////
	// Generate baseband signals
	////////////////////////////////////////////////////////////

	tstart = clock();

	// Update receiver time
	grx = incGpsTime(grx, 0.1);
    
    PeriodicFunctionTable<double, 1024> sin_table([](double t) -> double {return std::sin(t);}, 1.0);
    PeriodicFunctionTable<double, 1024> cos_table([](double t) -> double {return std::cos(t);}, 1.0);
    
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

	for (int iumd = 1; iumd < numd; iumd++)
	{
        // Per-channel gain, maximum value of 1.
        double gain[MAX_CHAN];

		for (int i = 0; i < MAX_CHAN; i++)
		{
			if (chan[i].prn>0)
			{
				// Refresh code phase and data bit counters
				range_t rho;
				const int sv = chan[i].prn-1;

				// Current pseudorange
				computeRange(&rho, eph[ieph][sv], &ionoutc, grx, xyz[0]);

				chan[i].azel[0] = rho.azel[0];
				chan[i].azel[1] = rho.azel[1];

				// Update code phase and data bit counters
				computeCodePhase(&chan[i], rho, 0.1);

				// Signal gain.
				gain[i] = constellation_gain.ComputeGain(
                    chan[i].prn, rho.d, rho.azel);
			}
		}

		// There is no intrinsic need to loop here, but we do so to avoid computing
		// the pseudorange too frequently.
		const int samples_in_100ms = std::lround((1.0 / 10.0) * samp_freq);
		for (int isamp = 0; isamp < samples_in_100ms; isamp++)
		{
			double i_acc = 0;
			double q_acc = 0;

			for (int i = 0; i < MAX_CHAN; i++)
			{
				if (chan[i].prn>0)
				{
                    const double coeff = chan[i].dataBit * chan[i].codeCA * gain[i];
                    const double ip = coeff * cos_table.LookupValue(chan[i].carr_phase);
                    const double qp = coeff * sin_table.LookupValue(chan[i].carr_phase);

					// Accumulate for all visible satellites
					i_acc += ip;
					q_acc += qp;

					// Update code phase
					chan[i].code_phase += chan[i].f_code * delt;

					if (chan[i].code_phase>=CA_SEQ_LEN)
					{
						chan[i].code_phase -= CA_SEQ_LEN;

						chan[i].icode++;
					
						if (chan[i].icode>=20) // 20 C/A codes = 1 navigation data bit
						{
							chan[i].icode = 0;
							chan[i].ibit++;
						
							if (chan[i].ibit>=30) // 30 navigation data bits = 1 word
							{
								chan[i].ibit = 0;
								chan[i].iword++;
                            }

							// Set new navigation data bit
							chan[i].dataBit = (int)((chan[i].dwrd[chan[i].iword]>>(29-chan[i].ibit)) & 0x1UL)*2-1;
						}
					}

					// Set current code chip
					chan[i].codeCA = chan[i].ca[(int)chan[i].code_phase]*2-1;

					// Update carrier phase
					chan[i].carr_phase += chan[i].f_carr * delt;

					if (chan[i].carr_phase >= 1.0)
						chan[i].carr_phase -= 1.0;
					else if (chan[i].carr_phase<0.0)
						chan[i].carr_phase += 1.0;
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
		// Update navigation message and channel allocation every 30 seconds
		//

        const int igrx = (int)(grx.sec*10.0+0.5);
		if (igrx % 300 == 0) // Every 30 seconds
		{
			// Update navigation message
			for (int i = 0; i < MAX_CHAN; i++)
			{
				if (chan[i].prn>0)
					generateNavMsg(grx, &chan[i], 0);
			}

			// Refresh ephemeris and subframes
			// Quick and dirty fix. Need more elegant way.
			for (int sv=0; sv < MAX_SAT; sv++)
			{
				if (eph[ieph+1][sv].vflg==1)
				{
					const double dt = subGpsTime(eph[ieph+1][sv].toc, grx);
					if (dt < SECONDS_IN_HOUR)
					{
						ieph++;
						for (int i = 0; i < MAX_CHAN; i++)
						{
							// Generate new subframes if allocated
							if (chan[i].prn!=0) 
								eph2sbf(eph[ieph][chan[i].prn-1], ionoutc, chan[i].sbf);
						}
					}
						
					break;
				}
			}

			// Update channel allocation
			allocateChannel(chan, eph[ieph], ionoutc, grx, xyz[0], elevation_mask_deg);

			// Show details about simulated channels
			if (verbose)
			{
				fprintf(stderr, "\n");
				for (int i = 0; i < MAX_CHAN; i++)
				{
					if (chan[i].prn>0)
						fprintf(stderr, "%02d %6.1f %5.1f %11.1f %5.1f\n", chan[i].prn,
							chan[i].azel[0]*R2D, chan[i].azel[1]*R2D, chan[i].rho0.d, chan[i].rho0.iono_delay);
				}
			}
		}

		// Update receiver time
		grx = incGpsTime(grx, 0.1);

		// Update time counter
		fprintf(stderr, "\rTime into run = %4.1f", subGpsTime(grx, g0));
		fflush(stdout);
	}

	tend = clock();

	fprintf(stderr, "\nDone!\n");

	// Process time
	fprintf(stderr, "Process time = %.1f [sec]\n", (double)(tend-tstart)/CLOCKS_PER_SEC);

	return(0);
}
