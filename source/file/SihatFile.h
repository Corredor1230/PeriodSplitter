#pragma once

#include<Windows.h>
#include<sndfile.h>
#include<commdlg.h>
#include<filesystem>
#include<shlobj.h>
#include"include/SitranoHeader.h"
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
		Analyzer ana(config);
		Sihat::Results r = ana.analyze(unit, settings);

		for (auto& amps : r.hResults.amps)
			if (!amps.empty()) Sihat::filterVector(amps, 4, true);
		for (auto& freqs : r.hResults.freqs)
			if (!freqs.empty()) Sihat::filterVector(freqs, 40, false);

		std::string csvName = Sihat::getRawFilename(filename);
		std::string dir = info.outDir.empty() ? "sihat" : info.outDir;

		std::string fName = info.prefix + "_" + csvName;
		std::string ext = "." + info.extension;

		if (settings.sourceSeparation) {
			Sihat::saveHarmonicDataSihat(r.hResults, r.stResults, config, r.pitch, dir, fName, ext);
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