#pragma once

#include<iostream>
#include<fstream>
#include<sndfile.h>
#include<vector>
#include<array>
#include<sndfile.h>
#include<string>

class Splitter
{
public:

	Splitter(SF_INFO& sf);
	~Splitter() {};

	void split(std::vector<float>& audioFile, std::vector<int>& splitSamples, std::string filename = "Audio",
		bool useNumOutputFiles = false, int numOutputFiles = 0);

	void writeCsvFile(std::vector<int>& outputInfo, std::string& filename);
	std::vector<std::pair<int, int>> loadCSV(const std::string& path, int totalSamples);

private:

	SF_INFO& sf;
	void writeAudioChunk(std::vector<float>& chunk, const std::string& filename);

	std::vector<int> clipList;

	int startSamp = 0;

};