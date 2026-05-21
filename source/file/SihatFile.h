#pragma once

#include<Windows.h>
#include<sndfile.h>
#include<commdlg.h>
#include<filesystem>
#include<shlobj.h>
#include<iostream>
#include"include/SitranoHeader.h"
#include"include/ResynthHeader.h"
#include"analysis/SitranoAnalysis.h"

namespace fs = std::filesystem;

namespace SihatFile {

	struct DirAndFiles {
		std::string directory;
		std::vector<std::string> files;
	};

	struct OutInfo {
		std::string outDir;
		std::string prefix;
		std::string extension;
	};

	struct AudioFile {
		int sampleRate = 44100;

	};

	inline std::string openFileDialog() {
		OPENFILENAMEA ofn;
		char szFile[260];
		ZeroMemory(&ofn, sizeof(ofn));
		ofn.lStructSize = sizeof(ofn);
		ofn.hwndOwner = nullptr;
		ofn.lpstrFile = szFile;
		ofn.lpstrFile[0] = '\0';
		ofn.nMaxFile = sizeof(szFile);
		ofn.lpstrFilter = "WAV Files\0*wav\0All Files\0*.*\0";
		ofn.nFilterIndex = 1;
		ofn.lpstrFileTitle = nullptr;
		ofn.nMaxFileTitle = 0;
		ofn.lpstrInitialDir = nullptr;
		ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

		if (GetOpenFileNameA(&ofn) == TRUE) {
			return std::string(ofn.lpstrFile);
		}

		return "";
	}

	inline std::string openJsonDialog() {
		OPENFILENAMEA ofn;
		char szFile[260];
		ZeroMemory(&ofn, sizeof(ofn));
		ofn.lStructSize = sizeof(ofn);
		ofn.hwndOwner = nullptr;
		ofn.lpstrFile = szFile;
		ofn.lpstrFile[0] = '\0';
		ofn.nMaxFile = sizeof(szFile);
		ofn.lpstrFilter = "JSON Files\0*json\0All Files\0*.*\0";
		ofn.nFilterIndex = 1;
		ofn.lpstrFileTitle = nullptr;
		ofn.nMaxFileTitle = 0;
		ofn.lpstrInitialDir = nullptr;
		ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

		if (GetOpenFileNameA(&ofn) == TRUE) {
			return std::string(ofn.lpstrFile);
		}

		return "";
	}

	inline std::string openSihatDialog() {
		OPENFILENAMEA ofn;
		char szFile[260];
		ZeroMemory(&ofn, sizeof(ofn));
		ofn.lStructSize = sizeof(ofn);
		ofn.hwndOwner = nullptr;
		ofn.lpstrFile = szFile;
		ofn.lpstrFile[0] = '\0';
		ofn.nMaxFile = sizeof(szFile);
		ofn.lpstrFilter = "sihat Files\0*sihat\0All Files\0*.*\0";
		ofn.nFilterIndex = 1;
		ofn.lpstrFileTitle = nullptr;
		ofn.nMaxFileTitle = 0;
		ofn.lpstrInitialDir = nullptr;
		ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

		if (GetOpenFileNameA(&ofn) == TRUE) {
			return std::string(ofn.lpstrFile);
		}

		return "";
	}

	/**
	 * Loads .sihat files and parses them into Sihat structs
	 */
	inline Synth::Sihat loadSihatFile(const std::string& fullpath)
	{
		std::ifstream inFile(fullpath, std::ios::binary);

		// 2. Always check if the file opened successfully
		if (!inFile.is_open()) {
			std::cerr << "Failed to open file for reading: " << fullpath << "\n";
			Synth::Sihat temp;
			return temp;
		}

		Synth::Sihat data;

		// Header load
		{
			inFile.read(reinterpret_cast<char*>(&data.header.sampleRate), sizeof(uint32_t));
			inFile.read(reinterpret_cast<char*>(&data.header.f0), sizeof(float));
			uint32_t stringSize = 0;
			inFile.read(reinterpret_cast<char*>(&stringSize), sizeof(uint32_t));
			if (stringSize > 0) {
				// Resize the string to allocate enough memory buffer space
				data.header.filename.resize(stringSize);

				// Read directly into the allocated memory buffer
				inFile.read(&data.header.filename[0], stringSize);
			}
		}

		//Block 1
		{
			inFile.read(reinterpret_cast<char*>(&data.harmonic.numFrames), sizeof(uint32_t));
			std::vector<uint32_t> tempIndex;
			tempIndex.resize(data.harmonic.numFrames);

			inFile.read(reinterpret_cast<char*>(tempIndex.data()), data.harmonic.numFrames * sizeof(uint32_t));

			data.harmonic.indices = tempIndex;
		}

		//Block 2: Amp
		{
			uint32_t numHarmonics = 0;
			inFile.read(reinterpret_cast<char*>(&numHarmonics), sizeof(uint32_t));
			data.harmonic.amp.resize(numHarmonics);
			for (int i = 0; i < numHarmonics; i++)
			{
				uint32_t numFrames = 0;
				inFile.read(reinterpret_cast<char*>(&numFrames), sizeof(uint32_t));
				std::vector<float> tempAmp;
				tempAmp.resize(numFrames);

				inFile.read(reinterpret_cast<char*>(tempAmp.data()), numFrames * sizeof(float));

				data.harmonic.amp[i] = tempAmp;
			}
		}

		// Block 3: Phase
		{
			uint32_t numHarmonics = 0;
			inFile.read(reinterpret_cast<char*>(&numHarmonics), sizeof(uint32_t));
			data.harmonic.pha.resize(numHarmonics);
			for (int i = 0; i < numHarmonics; i++)
			{
				uint32_t numFrames = 0;
				inFile.read(reinterpret_cast<char*>(&numFrames), sizeof(uint32_t));
				std::vector<float> tempPha;
				tempPha.resize(numFrames);

				inFile.read(reinterpret_cast<char*>(tempPha.data()), numFrames * sizeof(float));

				data.harmonic.pha[i] = tempPha;
			}
		}

		// Block 4: Frequency
		{
			uint32_t numHarmonics = 0;
			inFile.read(reinterpret_cast<char*>(&numHarmonics), sizeof(uint32_t));
			data.harmonic.freq.resize(numHarmonics);
			for (int i = 0; i < numHarmonics; i++)
			{
				uint32_t numFrames = 0;
				inFile.read(reinterpret_cast<char*>(&numFrames), sizeof(uint32_t));
				std::vector<float> tempFreq;
				tempFreq.resize(numFrames);

				inFile.read(reinterpret_cast<char*>(tempFreq.data()), numFrames * sizeof(float));

				data.harmonic.freq[i] = tempFreq;
			}
		}

		auto readFloatVector = [&inFile](std::vector<float>& vec) {
            uint32_t size = 0;
            inFile.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
            if (size > 0) {
                vec.resize(size);
                inFile.read(reinterpret_cast<char*>(vec.data()), size * sizeof(float));
            }
        };

        // Block 5
        {
            // 1. Read the Range
            inFile.read(reinterpret_cast<char*>(&data.transient.meta.tStart), sizeof(uint32_t));
            inFile.read(reinterpret_cast<char*>(&data.transient.meta.tEnd), sizeof(uint32_t));

            // Important parameters
            inFile.read(reinterpret_cast<char*>(&data.transient.meta.envHopSize), sizeof(uint32_t));
            inFile.read(reinterpret_cast<char*>(&data.transient.meta.specHopSize), sizeof(uint32_t));
            inFile.read(reinterpret_cast<char*>(&data.transient.meta.specWindowSize), sizeof(uint32_t));
            inFile.read(reinterpret_cast<char*>(&data.transient.meta.specNumBins), sizeof(uint32_t));

            // 2. Read scalar metrics
            inFile.read(reinterpret_cast<char*>(&data.transient.riseTime), sizeof(int32_t));
            inFile.read(reinterpret_cast<char*>(&data.transient.peakAmp), sizeof(float));

            // 3. Read 1D analysis envelopes
            uint32_t envSize = 0;
            inFile.read(reinterpret_cast<char*>(&envSize), sizeof(uint32_t));
			data.transient.envelope.env.resize(envSize);
            
            if (envSize > 0) {
                inFile.read(reinterpret_cast<char*>(&data.transient.envelope.firstIndex), sizeof(uint32_t));
                inFile.read(reinterpret_cast<char*>(data.transient.envelope.env.data()), envSize * sizeof(float));
            }

            readFloatVector(data.transient.centroid);
            readFloatVector(data.transient.flatness);

            // 4. Read the Band Partials
            inFile.read(reinterpret_cast<char*>(&data.transient.meta.numBands), sizeof(uint32_t));

            // 5. Read the flattened band envelopes
            readFloatVector(data.transient.bands);

            // Metadata for the harmonic section
            inFile.read(reinterpret_cast<char*>(&data.transient.meta.harmHopSize), sizeof(uint32_t));
            inFile.read(reinterpret_cast<char*>(&data.transient.meta.harmStartSample), sizeof(uint32_t));

            uint32_t numOvertones = 0;
            inFile.read(reinterpret_cast<char*>(&numOvertones), sizeof(uint32_t));
            data.transient.overtones.resize(numOvertones);

            // Harmonic section proper
            for (auto& over : data.transient.overtones) {
                // Assuming SpectralBin contains double freq, double mag, double amp
                inFile.read(reinterpret_cast<char*>(&over.target.freq), sizeof(double));
                inFile.read(reinterpret_cast<char*>(&over.target.mag), sizeof(double));
                inFile.read(reinterpret_cast<char*>(&over.target.amp), sizeof(double));

                uint32_t frameNum = 0;
                inFile.read(reinterpret_cast<char*>(&frameNum), sizeof(uint32_t));
                over.envelope.resize(frameNum);

                // Assuming EnvelopePoint contains float freq, float crestFactor
                for (auto& point : over.envelope) {
                    inFile.read(reinterpret_cast<char*>(&point.freq), sizeof(float));
                    inFile.read(reinterpret_cast<char*>(&point.crestFactor), sizeof(float));
                }
            }

            // Floor values
            readFloatVector(data.transient.floors);
        }

        // Block 5: Read the RMS values for each section
        {
            inFile.read(reinterpret_cast<char*>(&data.harmonic.rms), sizeof(float));
            inFile.read(reinterpret_cast<char*>(&data.transient.rms), sizeof(float));
        }

        return data; // Make sure to return your loaded struct!
	}

	inline void writeAudioFile(const std::vector<float>& audioData, const std::string& filename, const std::string& folderPath, uint32_t sampleRate = 96000) 
	{
		if (audioData.empty()) {
			std::cerr << "Error: Audio data is empty. Nothing to write.\n";
			return;
		}

		// Safely combine folder path and filename (C++17)
		std::filesystem::path fullPath = std::filesystem::path(folderPath) / filename;

		SF_INFO sfinfo;
		sfinfo.frames = audioData.size();
		sfinfo.samplerate = sampleRate;
		sfinfo.channels = 1; // Mono
		sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT; // 32-bit float WAV

		// Open the file for writing
		SNDFILE* outfile = sf_open(fullPath.string().c_str(), SFM_WRITE, &sfinfo);

		if (!outfile) {
			std::cerr << "Error opening file for writing: " << sf_strerror(nullptr) << "\n";
			return;
		}

		// Write the float data directly
		sf_count_t framesWritten = sf_write_float(outfile, audioData.data(), audioData.size());

		if (framesWritten != static_cast<sf_count_t>(audioData.size())) {
			std::cerr << "Warning: Wrote " << framesWritten << " frames out of " << audioData.size() << "\n";
		}

		sf_close(outfile);
	}

	/**
	 * Opens the Windows folder picker dialog (the same as the file picker).
	 * Returns the selected folder path, or an empty string if canceled.
	*/
	inline std::string openFolderDialog() {
		std::string folderPath = "";

		HRESULT hrInit = CoInitialize(NULL);
		if (FAILED(hrInit)) return folderPath;

		CoInitialize(NULL);

		IFileOpenDialog* pFileOpen;
		HRESULT hr = CoCreateInstance(__uuidof(FileOpenDialog), NULL, CLSCTX_ALL,
			IID_IFileOpenDialog, reinterpret_cast<void**>(&pFileOpen));

		if (SUCCEEDED(hr)) {

			DWORD dwOptions;
			hr = pFileOpen->GetOptions(&dwOptions);
			if (SUCCEEDED(hr)) {
				hr = pFileOpen->SetOptions(dwOptions | FOS_PICKFOLDERS);
			}

			hr = pFileOpen->Show(NULL);

			if (SUCCEEDED(hr)) {
				IShellItem* pItem;
				hr = pFileOpen->GetResult(&pItem);
				if (SUCCEEDED(hr)) {
					PWSTR pszFilePath;
					hr = pItem->GetDisplayName(SIGDN_FILESYSPATH, &pszFilePath);

					if (SUCCEEDED(hr)) {
						int bufferSize = WideCharToMultiByte(CP_UTF8, 0, pszFilePath, -1, NULL, 0, NULL, NULL);
						if (bufferSize > 0) {
							std::string mbPath(bufferSize, 0);
							WideCharToMultiByte(CP_UTF8, 0, pszFilePath, -1, &mbPath[0], bufferSize, NULL, NULL);
							mbPath.resize(bufferSize - 1);
							folderPath = mbPath;
						}
						CoTaskMemFree(pszFilePath);
					}
					pItem->Release();
				}
				pFileOpen->Release();
			}
			CoUninitialize();
			return folderPath;
		}
		return folderPath;
	}

	inline DirAndFiles getFileListFromExtension(const std::string& inDirectory,
		const std::string& extension) {
		std::vector<std::string> filenames;
		std::string directory = inDirectory;

		try {
			for (const auto& entry : fs::directory_iterator(directory)) {
				if (entry.is_regular_file() && entry.path().extension() == extension) {
					filenames.push_back(entry.path().filename().string());
				}
			}
		}
		catch (const fs::filesystem_error& e) {
			std::cerr << "Filesystem error: " << e.what() << '\n';
		}
		catch (const std::exception& e) {
			std::cerr << "Error: " << e.what() << '\n';
		}

		return DirAndFiles{ directory, filenames };
	}

	inline std::vector<float> getAudioFromFile(const std::string& filename,
		SF_INFO& sfInfo, float maxLength = 30.f){
		SNDFILE* file = sf_open(filename.c_str(), SFM_READ, &sfInfo);
		int maxNumSamples = sfInfo.samplerate * maxLength;
		const int numFrames = sfInfo.frames;
		int numSamps = sfInfo.frames > maxNumSamples ? maxNumSamples : sfInfo.frames;
		std::vector<std::vector<float>> delaced(sfInfo.channels, std::vector<float>(numSamps));
		std::vector<float> samples(sfInfo.frames * sfInfo.channels);
		sf_readf_float(file, samples.data(), sfInfo.frames);

		for (int chan = 0; chan < sfInfo.channels; chan++)
		{
			for (int samp = 0; samp < numSamps; samp++)
			{
				delaced[chan][samp] = samples[samp * sfInfo.channels + chan];
			}
		}

		return delaced[0];
	}

	/**
	* @brief processes an audio file with the selected config and settings
	* @param prefix is empty, but if included, it will be separated from the name by an underscore
	* @param extension Do not include the period, only the name of the extension. 
	*/
	inline void processFile(const std::string& filename,
		const OutInfo& info,
		const Sihat::AnalysisConfig& config,
		const Sihat::Settings& settings)
	{
		SF_INFO sfInfo;
		std::vector<float> delaced = SihatFile::getAudioFromFile(filename, sfInfo);
		Sihat::normalizeByMaxAbs(delaced);

		Sihat::AnalysisUnit unit{ delaced, filename, (float)sfInfo.samplerate };
		uint32_t sr = static_cast<uint32_t>(sfInfo.samplerate);
		Analyzer ana(config);
		Sihat::Results r = ana.analyze(unit, settings);

		for (auto& amps : r.hResults.amps)
			if (!amps.empty()) Sihat::filterVector(amps, 4, true);
		for (auto& freqs : r.hResults.freqs)
			if (!freqs.empty()) Sihat::filterVector(freqs, 40, false);

		std::string fName = Sihat::getRawFilename(filename);
		std::string prefix = info.prefix;
		std::string dir = info.outDir.empty() ? "sihat" : info.outDir;

		std::string ext = "." + info.extension;

		if (settings.sourceSeparation) {
			Sihat::saveHarmonicDataSihat(r.hResults, r.stResults, config, r.pitch, sr, dir, prefix, fName, ext);
		}
		else
		{
			Sihat::saveHarmonicDataSihat(r.hResults, r.tResults, r.pitch, dir, fName, ext);
		}
	}

	inline void processFolder(const std::string& inputDir,
		const OutInfo& i,
		const Sihat::AnalysisConfig& config,
		const Sihat::Settings& settings)
	{
		OutInfo info = i;
		SihatFile::DirAndFiles dnf = SihatFile::getFileListFromExtension(inputDir, ".wav");

		if (dnf.files.empty()) {
			std::cerr << "No .wav files found.\n";
			return;
		}

		for (const std::string& fn : dnf.files)
		{
			std::string filename = dnf.directory + "/" + fn;
			try {
				std::cout << "--- Processing: " << filename << " ---\n";
				processFile(filename, info, config, settings);
				std::cout << "--- Finished: " << filename << " ---\n";
			}
			catch (const std::exception& e) {
				std::cerr << "!!! FAILED " << filename << ": " << e.what() << "\n";
			}
		}

		std::cout << "Batch processing complete.\n";
	}

	inline void exportSeparatedAudio(const std::vector<float>& interleaved, const std::string& outputDir, int sampleRate = 96000, const std::string& name = "")
	{
		if (interleaved.size() % 2 != 0) {
			std::cerr << "Error: Interleaved audio vector size is not even. Cannot de-interleave." << std::endl;

			return;
		}

		size_t numFrames = interleaved.size() / 2;

		std::vector<float> hAudio(numFrames);
		std::vector<float> pAudio(numFrames);

		for (size_t i = 0; i < numFrames; i++){
			hAudio[i] = interleaved[i * 2];
			pAudio[i] = interleaved[i * 2 + 1];
		}

		SF_INFO sfinfo;
		sfinfo.channels = 1;
		sfinfo.samplerate = sampleRate;

		sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;

		std::filesystem::path dir(outputDir);
		if (!std::filesystem::exists(dir)){
			std::filesystem::create_directories(dir);
		}

		auto writeWav = [&](const std::vector<float>& audio, const std::string& filename) -> bool {
			std::filesystem::path filePath = dir / filename;

			SNDFILE* outfile = sf_open(filePath.string().c_str(), SFM_WRITE, &sfinfo);
			if (!outfile) {
				std::cerr << "Error: could not open" << filePath << " for writing." << std::endl;
				std::cerr << sf_strerror(NULL) << std::endl;
				return false;
			}

			sf_count_t framesWritten = sf_write_float(outfile, audio.data(), audio.size());
			if (framesWritten != audio.size()) {
				std::cerr << "Warning: Not all frames were written to " << filename << std::endl;
			}

			sf_close(outfile);
			return true;
		};

		bool successH = writeWav(hAudio, name + "harmonic.wav");
		bool successP = writeWav(pAudio, name + "percussive.wav");

		if (successH) std::cout << "Harmonic data successfully written." << std::endl;
		if (successP) std::cout << "Percussive data successfully written." << std::endl;
	}

	inline std::string getFilenameFromPath(const std::string& path)
	{
		return std::filesystem::path(path).stem().string();
	}
}