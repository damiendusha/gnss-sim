// Copyright (c) 2020 Damien Dusha
// SPDX-License-Identifier: MIT

#include "sample_writer.h"

#include <cmath>

namespace
{
double ClipSample(double sample)
{
    if (sample > 1.0)
        return 1.0;

    if (sample < -1.0)
        return -1.0;

    return sample;
}

int16_t ScaleSample(double sample)
{
    constexpr double kScale = 32767.0;  // 2^15 - 1;
    return static_cast<int16_t>(std::lround(sample * kScale));
}

}   // namespace

SampleWriter::~SampleWriter()
{
    FlushSamples();

    if (fp_ != nullptr)
    {
        fclose(fp_);
    }
}

bool SampleWriter::OpenFile(const std::string &filename) {
    if (fp_ != nullptr)
        return false;

    fp_ = fopen(filename.c_str(), "wb");
    return fp_ != nullptr;
}


bool SampleWriter::WriteSample(double i_sample, double q_sample)
{
    i_sample = ClipSample(i_sample);
    q_sample = ClipSample(q_sample);

    sample_buffer_[samples_stored_ * 2] = ScaleSample(i_sample);
    sample_buffer_[samples_stored_ * 2 + 1] = ScaleSample(q_sample);
    ++samples_stored_;

    return (samples_stored_ * 2) >= sample_buffer_.size() ?
        FlushSamples() : true;    
}

bool SampleWriter::FlushSamples()
{
    if (fp_ == nullptr)
    {
        samples_stored_ = 0;
        return false;
    }
    
    if (samples_stored_ == 0)
        return true;

    const int num_elements = 2 * samples_stored_;
    samples_stored_ = 0;
    return num_elements == fwrite(sample_buffer_.data(), sizeof(int16_t), num_elements, fp_);
}
