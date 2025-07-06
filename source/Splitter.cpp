#include"Splitter.h"

Splitter::Splitter(SF_INFO& sf) : sf(sf)
{

}

void Splitter::split(std::vector<float>& audioFile, std::vector<int>& splitSamples, std::string filename,
    bool useNumOutputFiles, int numOutputFiles)
{
    int outFiles = 0;
    if (useNumOutputFiles && numOutputFiles < splitSamples.size())
        outFiles = numOutputFiles;
    else
        outFiles = splitSamples.size() - 1;

    int lastSample = 0;
	for (int i = 0; i < outFiles; i++)
	{
        int beginSample = splitSamples[i];
        int endSample = splitSamples[i + 1];
        int chunkLength = endSample - beginSample;
        std::vector<float> chunk;

        if (chunkLength > 0)
        {
            chunk.resize(chunkLength);
            for (int samp = 0; samp < chunkLength; samp++)
            {
                chunk[samp] = audioFile[beginSample + samp];
            }
        }

        std::string outputFileName = filename + "_" + std::to_string(i) + ".wav";

		writeAudioChunk(chunk, outputFileName);
	}
}

void Splitter::writeCsvFile(std::vector<int>& outputInfo, std::string& filename)
{
    std::ofstream outputFile(filename, std::ios::out | std::ios::trunc);

    if (!outputFile.is_open())
    {
        std::cerr << "Error: could not open file" << filename << std::endl;
    }

    for (int i = 0; i < outputInfo.size(); i++)
    {
        std::string value = std::to_string(outputInfo[i]);
        outputFile << value << std::endl;
    }
}

void Splitter::writeAudioChunk(std::vector<float>& chunk, const std::string& filename)
{
    SF_INFO sfInfo = { 0 };
    sfInfo.samplerate = sf.samplerate;
    sfInfo.channels = sf.channels;
    sfInfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_24;

    SNDFILE* outFile = sf_open(filename.c_str(), SFM_WRITE, &sfInfo);
    if (!outFile) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    sf_write_float(outFile, chunk.data(), chunk.size());
    sf_close(outFile);

}