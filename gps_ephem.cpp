// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT
//
// Derived from https://github.com/osqzss/gps-sdr-sim (MIT Licence):
// Copyright (c) 2015-2020 Takuji Ebinuma

#include "gps_ephem.h"

#include <cmath>

void satpos(const ephem_t &eph, gpstime_t g, double *pos, double *vel, double *clk)
{
	// Computing Satellite Velocity using the Broadcast Ephemeris
	// http://www.ngs.noaa.gov/gps-toolbox/bc_velo.htm

	double tk = g.sec - eph.toe.sec;
	if(tk > SECONDS_IN_HALF_WEEK)
		tk -= SECONDS_IN_WEEK;
	else if(tk < -SECONDS_IN_HALF_WEEK)
		tk += SECONDS_IN_WEEK;

	const double mk = eph.m0 + eph.n*tk;
	double ek = mk;
	double ekold = ek + 1.0;
  
	double OneMinusecosE = 0; // Suppress the uninitialized warning.
	while (std::fabs(ek-ekold) > 1.0E-14)
	{
		ekold = ek;
		OneMinusecosE = 1.0-eph.ecc*cos(ekold);
		ek = ek + (mk-ekold+eph.ecc*sin(ekold))/OneMinusecosE;
	}

	const double sek = std::sin(ek);
	const double cek = std::cos(ek);

	const double ekdot = eph.n/OneMinusecosE;

	const double relativistic = -4.442807633E-10*eph.ecc*eph.sqrta*sek;

	const double pk = std::atan2(eph.sq1e2*sek,cek-eph.ecc) + eph.aop;
	const double pkdot = eph.sq1e2*ekdot/OneMinusecosE;

	const double s2pk = std::sin(2.0*pk);
	const double c2pk = std::cos(2.0*pk);

	const double uk = pk + eph.cus*s2pk + eph.cuc*c2pk;
	const double suk = std::sin(uk);
	const double cuk = std::cos(uk);
	const double ukdot = pkdot*(1.0 + 2.0*(eph.cus*c2pk - eph.cuc*s2pk));

	const double rk = eph.A*OneMinusecosE + eph.crc*c2pk + eph.crs*s2pk;
	const double rkdot = eph.A*eph.ecc*sek*ekdot + 2.0*pkdot*(eph.crs*c2pk - eph.crc*s2pk);

	const double ik = eph.inc0 + eph.idot*tk + eph.cic*c2pk + eph.cis*s2pk;
	const double sik = std::sin(ik);
	const double cik = std::cos(ik);
	const double ikdot = eph.idot + 2.0*pkdot*(eph.cis*c2pk - eph.cic*s2pk);

	const double xpk = rk*cuk;
	const double ypk = rk*suk;
	const double xpkdot = rkdot*cuk - ypk*ukdot;
	const double ypkdot = rkdot*suk + xpk*ukdot;

	const double ok = eph.omg0 + tk*eph.omgkdot - OMEGA_EARTH*eph.toe.sec;
	const double sok = std::sin(ok);
	const double cok = std::cos(ok);

	pos[0] = xpk*cok - ypk*cik*sok;
	pos[1] = xpk*sok + ypk*cik*cok;
	pos[2] = ypk*sik;

	const double tmp = ypkdot*cik - ypk*sik*ikdot;

	vel[0] = -eph.omgkdot*pos[1] + xpkdot*cok - tmp*sok;
	vel[1] = eph.omgkdot*pos[0] + xpkdot*sok + tmp*cok;
	vel[2] = ypk*cik*ikdot + ypkdot*sik;

	// Satellite clock correction
	tk = g.sec - eph.toc.sec;
	if(tk > SECONDS_IN_HALF_WEEK)
		tk -= SECONDS_IN_WEEK;
	else if(tk < -SECONDS_IN_HALF_WEEK)
		tk += SECONDS_IN_WEEK;

	clk[0] = eph.af0 + tk*(eph.af1 + tk*eph.af2) + relativistic - eph.tgd;  
	clk[1] = eph.af1 + 2.0*tk*eph.af2; 

	return;
}
