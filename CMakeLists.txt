# This file is part of:
# gnss-sim: A GNSS Signal Simulator
#
# Copyright (c) 2020 Damien Dusha
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(gnss_sim)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(gnss_sim 
    gpssim.cpp
    antenna_gain.cpp
    gps_channel.cpp
    gps_ephem.cpp
    gps_subframe.cpp
    gps_time.cpp
    geodesy.cpp
    ionosphere.cpp
    rinex2_reader.cpp
    sample_writer.cpp
    satellite_gain.cpp
)

target_link_libraries(gnss_sim m Eigen3::Eigen)

target_compile_options(gnss_sim PRIVATE -g -O3 -march=native -Wall -Wextra)

install(TARGETS gnss_sim DESTINATION bin)
