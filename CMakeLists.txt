# This file is part of:
# gnss-sim: A GNSS Signal Simulator
#
# Copyright (c) 2020 Damien Dusha
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.8)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(gnss_sim)

add_executable(gnss_sim 
    gpssim.cpp
    antenna_gain.cpp
    satellite_gain.cpp
    gps_ephem.cpp
    gps_subframe.cpp
    gps_time.cpp
    geodesy.cpp
    nmea_reader.cpp
    rinex2_reader.cpp
)

target_link_libraries(gnss_sim m)

target_compile_options(gnss_sim PRIVATE -O3 -march=native)

install(TARGETS gnss_sim DESTINATION bin)
