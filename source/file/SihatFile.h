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

	inline auto loadSihatHeader(std::ifstream& inFile, uint16_t vMajor, uint16_t vMinor) {
    // Automatically instantiates whatever type 'header' is inside the Sihat struct
    auto header = decltype(Synth::Sihat::header){}; 

		uint8_t tA, hA;
		inFile.read(reinterpret_cast<char*>(&tA), sizeof(uint8_t));
		header.type.transientAnalysis = (tA != 0);

		inFile.read(reinterpret_cast<char*>(&hA), sizeof(uint8_t));
		header.type.harmonicAnalysis = (hA != 0);

		inFile.read(reinterpret_cast<char*>(&header.sampleRate), sizeof(uint32_t));
		inFile.read(reinterpret_cast<char*>(&header.f0), sizeof(float));
		
		uint32_t stringSize = 0;
		inFile.read(reinterpret_cast<char*>(&stringSize), sizeof(uint32_t));
		if (stringSize > 0) {
			header.filename.resize(stringSize);
			inFile.read(&header.filename[0], stringSize);
		}
		
		return header;
	}

	inline auto loadSihatHarmonic(std::ifstream& inFile, uint16_t vMajor, uint16_t vMinor) {
		auto harmonic = decltype(Synth::Sihat::harmonic){};

		// Block 1: Indices
		inFile.read(reinterpret_cast<char*>(&harmonic.numFrames), sizeof(uint32_t));
		harmonic.indices.resize(harmonic.numFrames);
		inFile.read(reinterpret_cast<char*>(harmonic.indices.data()), harmonic.numFrames * sizeof(uint32_t));

		// Block 2: Amp
		uint32_t numHarmonics = 0;
		inFile.read(reinterpret_cast<char*>(&numHarmonics), sizeof(uint32_t));
		harmonic.amp.resize(numHarmonics);
		for (int i = 0; i < numHarmonics; i++) {
			uint32_t numFrames = 0;
			inFile.read(reinterpret_cast<char*>(&numFrames), sizeof(uint32_t));
			harmonic.amp[i].resize(numFrames);
			inFile.read(reinterpret_cast<char*>(harmonic.amp[i].data()), numFrames * sizeof(float));
		}

		// Block 3: Phase
		inFile.read(reinterpret_cast<char*>(&numHarmonics), sizeof(uint32_t));
		harmonic.pha.resize(numHarmonics);
		for (int i = 0; i < numHarmonics; i++) {
			uint32_t numFrames = 0;
			inFile.read(reinterpret_cast<char*>(&numFrames), sizeof(uint32_t));
			harmonic.pha[i].resize(numFrames);
			inFile.read(reinterpret_cast<char*>(harmonic.pha[i].data()), numFrames * sizeof(float));
		}

		// Block 4: Frequency
		inFile.read(reinterpret_cast<char*>(&numHarmonics), sizeof(uint32_t));
		harmonic.fRatio.resize(numHarmonics);
		for (int i = 0; i < numHarmonics; i++) {
			uint32_t numFrames = 0;
			inFile.read(reinterpret_cast<char*>(&numFrames), sizeof(uint32_t));
			harmonic.fRatio[i].resize(numFrames);
			inFile.read(reinterpret_cast<char*>(harmonic.fRatio[i].data()), numFrames * sizeof(float));
		}

		inFile.read(reinterpret_cast<char*>(&harmonic.rms), sizeof(float));

		return harmonic;
	}

	inline auto loadSihatTransient(std::ifstream& inFile, uint16_t vMajor, uint16_t vMinor) {
		auto transient = decltype(Synth::Sihat::transient){};

		auto readFloatVector = [&inFile](std::vector<float>& vec) {
			uint32_t size = 0;
			inFile.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
			if (size > 0) {
				vec.resize(size);
				inFile.read(reinterpret_cast<char*>(vec.data()), size * sizeof(float));
			}
		};

		// 1. Read the Range
		inFile.read(reinterpret_cast<char*>(&transient.meta.tStart), sizeof(uint32_t));
		inFile.read(reinterpret_cast<char*>(&transient.meta.tEnd), sizeof(uint32_t));

		// Important parameters
		inFile.read(reinterpret_cast<char*>(&transient.meta.envHopSize), sizeof(uint32_t));
		inFile.read(reinterpret_cast<char*>(&transient.meta.specHopSize), sizeof(uint32_t));
		inFile.read(reinterpret_cast<char*>(&transient.meta.floorHopSize), sizeof(uint32_t));
		inFile.read(reinterpret_cast<char*>(&transient.meta.specWindowSize), sizeof(uint32_t));
		inFile.read(reinterpret_cast<char*>(&transient.meta.specNumBins), sizeof(uint32_t));
		inFile.read(reinterpret_cast<char*>(&transient.meta.specNumFrames), sizeof(uint32_t));
		inFile.read(reinterpret_cast<char*>(&transient.meta.specNfft), sizeof(uint32_t));

		// 2. Read scalar metrics
		inFile.read(reinterpret_cast<char*>(&transient.riseTime), sizeof(int32_t));
		inFile.read(reinterpret_cast<char*>(&transient.peakAmp), sizeof(float));

		// 3. Read 1D analysis envelopes
		uint32_t envSize = 0;
		inFile.read(reinterpret_cast<char*>(&envSize), sizeof(uint32_t));
		transient.envelope.env.resize(envSize);
		
		if (envSize > 0) {
			inFile.read(reinterpret_cast<char*>(&transient.envelope.firstIndex), sizeof(uint32_t));
			inFile.read(reinterpret_cast<char*>(transient.envelope.env.data()), envSize * sizeof(float));
		}

		readFloatVector(transient.centroid);
		readFloatVector(transient.flatness);

		// 4. Read the Band Partials
		inFile.read(reinterpret_cast<char*>(&transient.meta.numBands), sizeof(uint32_t));

		// 5. Read the flattened band envelopes
		readFloatVector(transient.bands);

		// 6. Read the modes
		inFile.read(reinterpret_cast<char*>(&transient.tModes.startInd), sizeof(uint32_t));
		inFile.read(reinterpret_cast<char*>(&transient.tModes.length), sizeof(uint32_t));
		uint32_t modeNum = 0;
		inFile.read(reinterpret_cast<char*>(&modeNum), sizeof(uint32_t));
		transient.tModes.modes.resize(modeNum);
		
		for (auto& mode : transient.tModes.modes) {
			inFile.read(reinterpret_cast<char*>(&mode.freq), sizeof(float));
			inFile.read(reinterpret_cast<char*>(&mode.amp), sizeof(float));
			inFile.read(reinterpret_cast<char*>(&mode.phase), sizeof(float));
			inFile.read(reinterpret_cast<char*>(&mode.decay), sizeof(float));
		}

		// Metadata for the harmonic section
		inFile.read(reinterpret_cast<char*>(&transient.meta.harmHopSize), sizeof(uint32_t));
		inFile.read(reinterpret_cast<char*>(&transient.meta.harmStartSample), sizeof(uint32_t));

		uint32_t numOvertones = 0;
		inFile.read(reinterpret_cast<char*>(&numOvertones), sizeof(uint32_t));
		transient.overtones.resize(numOvertones);

		// Harmonic section proper
		for (auto& over : transient.overtones) {
			inFile.read(reinterpret_cast<char*>(&over.target.freq), sizeof(double));
			inFile.read(reinterpret_cast<char*>(&over.target.mag), sizeof(double));
			inFile.read(reinterpret_cast<char*>(&over.target.amp), sizeof(double));
			inFile.read(reinterpret_cast<char*>(&over.target.phase), sizeof(double));

			uint32_t frameNum = 0;
			inFile.read(reinterpret_cast<char*>(&frameNum), sizeof(uint32_t));
			over.envelope.resize(frameNum);

			for (auto& point : over.envelope) {
				inFile.read(reinterpret_cast<char*>(&point.freq), sizeof(float));
				inFile.read(reinterpret_cast<char*>(&point.crestFactor), sizeof(float));
				inFile.read(reinterpret_cast<char*>(&point.amp), sizeof(float));
			}
		}

		// Floor values
		readFloatVector(transient.floors);

		inFile.read(reinterpret_cast<char*>(&transient.rms), sizeof(float));

		return transient;
	}

	/**
	 * Loads .sihat files and parses them into Sihat structs
	 */
	inline Synth::Sihat loadSihatFile(const std::string& fullpath)
	{
		std::ifstream inFile(fullpath, std::ios::binary);

		if (!inFile.is_open()) {
			std::cerr << "Failed to open file for reading: " << fullpath << "\n";
			return Synth::Sihat{};
		}

		// 1. Verify Magic String
		char magic[4];
		inFile.read(magic, 4);
		if (inFile.gcount() != 4 || std::strncmp(magic, Synth::SIHAT_MAGIC, 4) != 0) {
			std::cerr << "[Error] File is not a valid Sihat format.\n";
			return Synth::Sihat{};
		}

		// 2. Read File Version
		uint16_t fileMajor = 0;
		uint16_t fileMinor = 0;
		inFile.read(reinterpret_cast<char*>(&fileMajor), sizeof(uint16_t));
		inFile.read(reinterpret_cast<char*>(&fileMinor), sizeof(uint16_t));

		// 3. Version Compatibility Check
		if (fileMajor > Synth::SIHAT_VERSION_MAJOR || 
		(fileMajor == Synth::SIHAT_VERSION_MAJOR && fileMinor > Synth::SIHAT_VERSION_MINOR)) {
			
			std::cerr << "[Error] File version v" << fileMajor << "." << fileMinor 
					<< " is newer than supported version v" 
					<< Synth::SIHAT_VERSION_MAJOR << "." << Synth::SIHAT_VERSION_MINOR 
					<< ". Please update your software.\n";
			return Synth::Sihat{};
		}

		// 4. Load Data
		Synth::Sihat data;
		
		// Notice how clean the data assignment is now
		data.header = loadSihatHeader(inFile, fileMajor, fileMinor);

		if (data.header.type.harmonicAnalysis) {
			data.harmonic = loadSihatHarmonic(inFile, fileMajor, fileMinor);
		}

		if (data.header.type.transientAnalysis) {
			data.transient = loadSihatTransient(inFile, fileMajor, fileMinor);
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
			Sihat::saveHarmonicDataSihat(r.hResults, r.stResults, config, settings, r.pitch, sr, dir, prefix, fName, ext);
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