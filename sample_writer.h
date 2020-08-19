// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT

#ifndef SAMPLE_BUFFER_H__
#define SAMPLE_BUFFER_H__

#include <string>
#include <cstdio>
#include <vector>
#include <cstdint>

class SampleWriter
{
  public:
    /// \brief SampleWriter Constructor
    ///
    /// \arg buffer_size_ is the number of 16-bit complex samples to write
    ///      in each chunk.
    SampleWriter(int buffer_size = 1024*1024)
    {
        sample_buffer_.resize(buffer_size * 2);
    }
        
    ~SampleWriter();

    bool OpenFile(const std::string &filename);

    // Writes a sample to the output file. The values are expected to be in the
    // range of -1.0 to +1.0. Values beyond these are clipped. These values
    // will be mapped across the full range of a signed 16-bit integer.
    bool WriteSample(double i_sample, double q_sample);

  private: 
    FILE* fp_ = nullptr;

    int samples_stored_ = 0;
    std::vector<int16_t> sample_buffer_;

    // Flushed all samples from the buffer to the file.
    bool FlushSamples();
};


#endif  // SAMPLE_BUFFER_H__

