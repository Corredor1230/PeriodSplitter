#include "NoiseTracker.h"

NoiseTracker::NoiseTracker(std::vector<float> audio, float start, float sampleRate, int fftSize = 8192, int hopSize = 64):
	sf(audio),
	startFreq(start),
	sr(sampleRate),
	N(fftSize),
	hop(hopSize)
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

void NoiseTracker::analyze() {
	int numFrames = useFrameTable ? frameTable.size() : (sf.size() - N) / hop;
}