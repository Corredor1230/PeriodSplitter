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

void Splitter::writeCsvFile(std::vector<std::vector<float>>& outInfo, std::string& filename,
    std::vector<float>& headers)
{
    std::ofstream outputFile(filename, std::ios::out | std::ios::trunc);

    if (!outputFile.is_open())
    {
        std::cerr << "Error: could not open file" << filename << std::endl;
    }

    for (int i = 0; i < outInfo.size(); i++)
    {
        outputFile << headers[i];
        if (i < headers.size() - 1) outputFile << ",";
    }
    outputFile << '\n';

    size_t numRows = outInfo[0].size();

    for (size_t row = 0; row < numRows; row++) {
        for (size_t col = 0; col < outInfo.size(); col++) {
            outputFile << outInfo[col][row];
            if (col < outInfo.size() - 1) outputFile << ",";
        }
        outputFile << '\n';
    }

    outputFile.close();
}

void writeCsvFile(std::vector<std::vector<float>>& outInfo, std::string& filename,
    std::vector<float>& columnHeaders, std::vector<int>& rowHeaders)
{
    std::ofstream out(filename);

    if (out.is_open())
    {
        std::cerr << "Error" << filename << std::endl;
    }

    out << " ,";

    for (int i = 0; i < columnHeaders.size(); i++)
    {
        out << columnHeaders[i];
        if (i + 1 < columnHeaders.size()) out << " ,";
    }

    out << '\n';

    for (size_t i = 0; i < outInfo.size(); ++i) {
        if (i < rowHeaders.size()) {
            out << rowHeaders[i] << ",";
        }
        else {
            out << " ,";  // safeguard: empty if no row header
        }

        for (size_t j = 0; j < outInfo[i].size(); ++j) {
            out << std::fixed << std::setprecision(3) << outInfo[i][j];
            if (j + 1 < outInfo[i].size()) out << ",";
        }
        out << "\n";
    }

    out.close();
}

std::vector<std::pair<int, int>> Splitter::loadCSV(const std::string& path, int totalSamples) 
{
    std::ifstream file(path);
    std::vector<int> indices;
    std::string line;

    while (std::getline(file, line)) {
        int index;
        if (sscanf(line.c_str(), "%d", &index) == 1)
            indices.push_back(index);
    }

    std::vector<std::pair<int, int>> result;
    for (size_t i = 0; i < indices.size(); ++i) {
        int start = indices[i];
        int end = (i + 1 < indices.size()) ? (indices[i + 1] - 1) : (totalSamples - 1);
        result.emplace_back(start, end);
    }

    return result;
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