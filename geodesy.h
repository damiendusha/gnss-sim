// This file is part of:
// gnss-sim: A GNSS Signal Simulator
//
// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT
//
// Derived from https://github.com/osqzss/gps-sdr-sim (MIT Licence):
// Copyright (c) 2015-2020 Takuji Ebinuma

#pragma once

#include <cmath>

class AzimuthElevation {
  public:
    AzimuthElevation() = default;

    /// \brief Creates an AzimuthElevation from degrees.
    static AzimuthElevation FromDegrees(double az_deg, double el_deg) {
        return AzimuthElevation(az_deg * (M_PI / 180.0), el_deg * (M_PI / 180.0));
    }

    /// \brief Creates an AzimuthElevation from radians.
    static AzimuthElevation FromRadians(double az_rad, double el_rad) {
        return AzimuthElevation(az_rad, el_rad);
    }

    void set_azimuth_deg(double az_deg) { set_azimuth_rad(az_deg * (M_PI / 180.0)); }
    void set_elevation_deg(double el_deg) { set_azimuth_rad(el_deg * (M_PI / 180.0)); }

    double azimuth_deg() const { return az_rad_ * (180.0 / M_PI); }
    double elevation_deg() const { return el_rad_ * (180.0 / M_PI); }
    
    void set_azimuth_rad(double az_rad) { az_rad_ = az_rad; }
    void set_elevation_rad(double el_rad) { el_rad_ = el_rad; }

    double azimuth_rad() const { return az_rad_; }
    double elevation_rad() const { return el_rad_; }

  private:
    constexpr AzimuthElevation(double az_rad, double el_rad) 
        : az_rad_(az_rad)
        , el_rad_(el_rad)
    {}

    double az_rad_ = 0;
    double el_rad_ = 0;
};

class GeodeticPosition {
  public:
    GeodeticPosition() = default;

    static GeodeticPosition FromDegrees(double lat_deg, double lon_deg, double height_m) {
        return GeodeticPosition(lat_deg * (M_PI / 180.0),
                                lon_deg * (M_PI / 180.0),
                                height_m);
    }

    static GeodeticPosition FromRadians(double lat_rad, double lon_rad, double height_m) {
        return GeodeticPosition(lat_rad, lon_rad, height_m);
    }

    void set_latitude_deg(double lat_deg) { set_latitude_rad(lat_deg * (M_PI / 180.0)); }
    void set_longitude_deg(double lon_deg) { set_longitude_rad(lon_deg * (M_PI / 180.0)); }
    void set_height_m(double height_m) { height_m_ = height_m; }

    double latitude_deg() const { return lat_rad_ * (180.0 / M_PI); }
    double longitude_deg() const { return lat_rad_ * (180.0 / M_PI); }
    double height_m() const { return height_m_; }
    
    void set_latitude_rad(double lat_rad) { lat_rad_ = lat_rad; }
    void set_longitude_rad(double lon_rad) { lon_rad_ = lon_rad; }

    double latitude_rad() const { return lat_rad_; }
    double longitude_rad() const { return lon_rad_; }    

  private:
    double lat_rad_ = 0;
    double lon_rad_ = 0;
    double height_m_ = 0;
    
    GeodeticPosition(double lat_rad, double lon_rad, double height_m)
        : lat_rad_(lat_rad)
        , lon_rad_(lon_rad)
        , height_m_(height_m)
    {
    }
};


/*! \brief Convert Earth-centered Earth-fixed (ECEF) into Lat/Lon/Height
 *  \param[in] xyz Input Array of X, Y and Z ECEF coordinates
 *  \param[out] llh Output Array of Latitude, Longitude and Height
 */
GeodeticPosition xyz2llh(const double *xyz);

/*! \brief Convert Lat/Long/Height into Earth-centered Earth-fixed (ECEF)
 *  \param[in] llh Input Array of Latitude, Longitude and Height
 *  \param[out] xyz Output Array of X, Y and Z ECEF coordinates
 */
void llh2xyz(const GeodeticPosition &llh, double *xyz);

/*! \brief Convert North-Eeast-Up to Azimuth + Elevation
 *  \param[in] neu Input position in North-East-Up format
 *  \returns Azimuth and elevation.
 */
AzimuthElevation neu2azel(const double *neu);

/*! \brief Convert Earth-centered Earth-Fixed to ?
 *  \param[in] xyz Input position as vector in ECEF format
 *  \param[in] t Intermediate matrix computed by \ref ltcmat
 *  \param[out] neu Output position as North-East-Up format
 */
void ecef2neu(const double *xyz, double t[3][3], double *neu);

/*! \brief Compute the intermediate matrix for LLH to ECEF
 *  \param[in] llh Input position in Latitude-Longitude-Height format
 *  \param[out] t Three-by-Three output matrix
 */
void ltcmat(const GeodeticPosition &llh, double t[3][3]);
