#pragma once

#include<Windows.h>
#include<sndfile.h>
#include<commdlg.h>
#include<filesystem>
#include<shlobj.h>
#include"include/json.hpp"

namespace fs = std::filesystem;

namespace SihatFile {

	struct DirAndFiles {
		std::string directory;
		std::vector<std::string> files;
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

	/**
	 * Opens the Windows folder picker dialog (the same as the file picker).
	 * Returns the selected folder path, or an empty string if canceled.
	*/
	inline std::string openFolderDialog() {
		std::string folderPath = "";
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

}