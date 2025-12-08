# SIHAT: Sitrano-Inspired Harmonic Analysis Tool

[![C++](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://isocpp.org/)
[![Build](https://img.shields.io/badge/build-passing-brightgreen.svg)](https://github.com/Corredor1230/PeriodSplitter)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

SIHAT is a high-performance C++ library and command-line tool for in-depth analysis of musical audio. It specializes in deconstructing a sound into its constituent partials, tracking the time-varying amplitude and frequency of each one.

This tool is ideal for academic research in musical acoustics, development of synthesis models (e.g., spectral or physical modeling), or any application requiring detailed, time-varying harmonic data. 

It can output highly detailed and low-dimensional data while retaining the distinctive harmonic features of single-note audio samples.

## Features

* **Pitch Detection:** Robust fundamental frequency ($f_0$) estimation.
* **Transient Detection:** Uses a sliding RMS window to detect the initial onset of the signal, and an autocorrelation algorithm to detect the end of the transient as the point where signal periodicity begins. 
* **Period Segmentation:** Intelligently segments the audio by pitch periods from the end of the initial transient to the point where all overtones' amplitudes fall below the specified threshold in decibels.
* **Overtone Analysis:** Identifies the most representative partials in the signal using a long FFT analysis.
* **Harmonic Tracking:** Tracks the precise frequency and amplitude of each partial over the entire duration of the audio file.
* **Highly Configurable:** Provides detailed control over every stage of the analysis pipeline.
* **Binary Output:** Saves analysis data in binary formats for easy loading into other environments (e.g., Python, MATLAB, or back into C++).

---

## Analysis Pipeline

The `Analyzer` class processes an audio file in a sequential, multi-stage pipeline. Each stage's output feeds into the next, allowing for a highly detailed and accurate analysis. Steps can also be turned off in which case their results will be replaced by default behaviors. Although this will yield less-detailed results, it will speed up performance significantly. 

1.  **Pitch Analysis (`PitchFinder`)**
    * The `PitchFinder` first analyzes the signal to determine a single, representative fundamental frequency ($f_0$) for the entire file.
    * This $f_0$ is used to guide the subsequent period and harmonic analysis stages.

2.  **Period Analysis (`PeriodCutter`)**
    * Instead of using a standard, fixed-size STFT, SIHAT segments the audio based on pitch periods.
    * It uses the `CorrelationSettings` (e.g., `transientRms`, `tFactor`) to detect and isolate transients, focusing the analysis on the stable, periodic portions of the sound.
    * The output is a list of sample indices (`sampleList`) marking the beginning of each valid analysis frame.
    * Standard FFT analysis is also possible using a different FFT window strategy.

3.  **Overtone Analysis (`OvertoneFinder`)**
    * This stage determines *what* to track. It analyzes an initial segment of the audio (defined by `overtoneStart`) to identify the most prominent overtones (partials).
    * This is crucial for inharmonic sounds (like pianos or guitars), as it finds the *actual* frequency of the partials rather than assuming they are perfect integer multiples of the fundamental.

4.  **Harmonic Tracking (`HarmonicTracker`)**
    * This is the final and most intensive stage. The `HarmonicTracker` takes the list of analysis frames (`sampleList`) and the list of partials to find (`topFreqs`).
    * For each frame, it performs a high-resolution spectral analysis to find the precise amplitude and frequency of each partial.
    * The result is a time-series dataset containing the amplitude and frequency "envelopes" for every partial.

---

## Integrated GUI

Using the [Dear ImGui](https://github.com/ocornut/imgui) library, you can edit the analysis parameters flexibly without having to alter the application's code. 

SIHAT's code follows the model-view-controller (MVC) paradigm, which means you can easily replace the GUI library for other common alternatives such as JUCE or QT and have the same DSP performance. The main DSP functions can be found within the EditorViewModel.h file. These are the only actions that need to be triggered from the GUI directly.

---

## Building from Source

This project is built using CMake.

### Dependencies

* [**CMake**](https://cmake.org/) (version 3.15 or higher)
* [**libsndfile**](http://www.mega-nerd.com/libsndfile/): For loading audio files.
* [**FFTW3**](https://www.fftw.org/): For high-performance Fast Fourier Transforms.
* [**PYin**](https://code.soundsoftware.ac.uk/projects/pyin): For pitch recognition.
* [**Dear ImGui**](https://github.com/ocornut/imgui): For the GUI

### Build Steps

1.  Clone the repository:
    ```bash
    git clone [https://github.com/Corredor1230/PeriodSplitter.git](https://github.com/Corredor1230/PeriodSplitter.git)
    cd SIHAT
    ```

2.  Create and navigate to the build directory:
    ```bash
    mkdir build
    cd build
    ```

3.  Run CMake and build the project:
    ```bash
    cmake ..
    make
    ```
    (or `cmake --build .` for a more general command)

---

## Configuration

Analysis parameters are set in `main.cpp` within the `Sitrano::AnalysisConfig` struct.

| Parameter | Struct | Description |
| :--- | :--- | :--- |
| `maxHarmonics` | `AnalysisConfig` | The maximum number of partials to search for and track. |
| `N` | `AnalysisConfig` | The default FFT size for spectral analysis. |
| `hopSize` | `AnalysisConfig` | The hop size to use if period analysis is disabled. |
| `tolerance` | `AnalysisConfig` | A general frequency tolerance (in Hz) for peak matching. |
| `verbose` | `AnalysisConfig` | Print detailed analysis progress to the console. |
| | | |
| `modeThreshold` | `PitchSettings` | Threshold for finding the mode of the pitch array. |
| `toleranceInCents` | `PitchSettings` | Tolerance value in cents to allow for differences between the detected pitch and the metadata embedded pitch value.
| `minFreq` / `maxFreq` | `PitchSettings` | The valid frequency range for the initial pitch search. |
| | | |
| `transientRms` | `CorrelationSettings` | The window size (in ms) for RMS-based transient detection. |
| `tFactor` | `CorrelationSettings` | The multiplier for transient detection (e.g., a 3.0x increase in RMS is a transient). |
| `periodOffset` | `CorrelationSettings` | How far (in ms) after a transient to begin period analysis. |
| `corrThreshold` | `CorrelationSettings` | The auto-correlation threshold (0.0-1.0) required for a frame to be considered "periodic." |
| | | |
| `overtoneStart` | `OvertoneSettings` | The sample index to start the initial overtone analysis from. |
| `fftSize` | `OvertoneSettings` | A specific (often larger) FFT size for the high-resolution overtone search. |
| | | |
| `windowStyle` | `HarmonicSettings` | The windowing method to use for harmonic tracking (e.g., `audioChunk`). |
| `oTolerance` | `HarmonicSettings` | A specific frequency tolerance (in Hz) for the final harmonic tracking stage. |

---

## Usage

Currently, the analysis is configured from the GUI, and then compiled, and run from the `SihatApplication()` class.

1.  **Compile and run:**
    ```bash
    cd build
    cmake ..
    ```

2.  **Configure your analysis:**
    Edit the values and variables in the GUI to produce .sihat files that contain the main information extracted from your audio file. You can process single audio files or batch process entire directories.

---

## Output

The program saves the analysis data to one binary file with extension .sihat:

* **`GEN_FILENAME.sihat`**: 

1. F0 (fundamental frequency) `float` containing the floating point value found for the input signal.

2. A `std::vector<int>` containing the sample index for the start of each analysis frame. The size of this vector defines the time dimension of the other files.

3. A `std::vector<std::vector<float>>` (Harmonic x Time) containing the tracked amplitude for each partial in each frame.

4. A `std::vector<std::vector<float>>` (Harmonic x Time) containing the tracked frequency for each partial in each frame.

5. A `std::vector<std::vector<float>>` (Harmonic x Time) containing the tracked ratio of overtone frequency to fundamental frequency ($f_0$) for each partial in each frame.

This file can be loaded by other C++ applications or read in scripting languages like Python (using `numpy.fromfile`) for plotting and further modeling.

## License

This project is licensed under the GPL v3 License - see the [LICENSE.md](LICENSE.md) file for details.