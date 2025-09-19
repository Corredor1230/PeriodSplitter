#include "NoiseTracker.h"

NoiseTracker::NoiseTracker(float start, float sampleRate):
	startFreq(start),
	sr(sampleRate)
{
	float f0 = startFreq;
	while (f0 < sr / 2.0) {
		float f_low = f0 / pow(2.0, 1.0 / 6.0);
		float f_high = f0 * pow(2.0, 1.0 / 6.0);
		bands.push_back({ f_low, f_high, f0 });
		f0 *= pow(2.0, 1.0 / 3.0);
	}
}

void NoiseTracker::applyFrameTable(std::vector<int> table)
{
	frameTable = table;
	useFrameTable = true;
}

void NoiseTracker::analyze(fftw_complex* out) {
	int numFrames = useFrameTable ? frameTable.size() : 
	results.assign()
}