#include"WTransform.h"

WTransform::WTransform(int windowSize)
: wSize(windowSize)
{
    plan, input, output = nullptr;

    initFFTW();
}

WTransform::~WTransform()
{
    if (plan != nullptr) fftwf_destroy_plan(plan);
    if (input != nullptr) fftwf_free(input);
    if (output != nullptr) fftwf_free(output);
}

void WTransform::initFFTW()
{
    input = (float*)fftwf_malloc(sizeof(float) * wSize);
    output = (fftwf_complex*)fftwf_malloc(sizeof(float) * wSize / 2 + 1);
    plan = fftwf_plan_dft_r2c_1d(wSize, input, output, FFTW_MEASURE);
}

std::vector<float> WTransform::wave(std::vector<float>& input, int windowSize)
{

}