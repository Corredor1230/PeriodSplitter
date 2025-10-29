# SIHAT: Sitrano-Inspired Harmonic Analysis Tool

[![C++](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://isocpp.org/)
[![Build](https://img.shields.io/badge/build-passing-brightgreen.svg)](https://github.com/YOUR_USERNAME/YOUR_REPO)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

SIHAT is a high-performance C++ library and command-line tool for in-depth analysis of musical audio. It specializes in deconstructing a sound into its constituent partials, tracking the time-varying amplitude and frequency of each one.

This tool is ideal for academic research in musical acoustics, development of synthesis models (e.g., spectral or physical modeling), or any application requiring detailed, time-varying harmonic data.

## Features

* **Pitch Detection:** Robust fundamental frequency ($f_0$) estimation.
* **Period Segmentation:** Intelligently segments the audio by pitch periods, with transient detection to find stable, periodic sections.
* **Inharmonic Overtone Analysis:** Identifies the most prominent partials in the signal, which are not assumed to be perfectly harmonic.
* **Harmonic Tracking:** Tracks the precise frequency and amplitude of each partial over the entire duration of the audio file.
* **Highly Configurable:** Provides detailed control over every stage of the analysis pipeline.
* **Binary Output:** Saves analysis data in efficient binary formats for easy loading into other environments (e.g., Python, MATLAB, or back into C++).

---

## Analysis Pipeline

The `Analyzer` class processes an audio file in a sequential, multi-stage pipeline. Each stage's output feeds into the next, allowing for a highly detailed and accurate analysis.

1.  **Pitch Analysis (`PitchFinder`)**
    * The `PitchFinder` first analyzes the signal to determine a single, representative fundamental frequency ($f_0$) for the entire file.
    * This $f_0$ is used to guide the subsequent period and harmonic analysis stages.

2.  **Period Analysis (`PeriodCutter`)**
    * Instead of using a standard, fixed-size STFT, Sitrano segments the audio based on pitch periods.
    * It uses the `CorrelationSettings` (e.g., `transientRms`, `tFactor`) to detect and isolate transients, focusing the analysis on the stable, periodic portions of the sound.
    * The output is a list of sample indices (`sampleList`) marking the beginning of each valid analysis frame.

3.  **Overtone Analysis (`OvertoneFinder`)**
    * This stage determines *what* to track. It analyzes an initial segment of the audio (defined by `overtoneStart`) to identify the most prominent overtones (partials).
    * This is crucial for inharmonic sounds (like pianos or guitars), as it finds the *actual* frequency of the partials rather than assuming they are perfect integer multiples of the fundamental.

4.  **Harmonic Tracking (`HarmonicTracker`)**
    * This is the final and most intensive stage. The `HarmonicTracker` takes the list of analysis frames (`sampleList`) and the list of partials to find (`topFreqs`).
    * For each frame, it performs a high-resolution spectral analysis to find the precise amplitude and frequency of each partial.
    * The result is a time-series dataset containing the amplitude and frequency "envelopes" for every partial.

---

## Building from Source

This project is built using CMake.

### Dependencies

* [**CMake**](https://cmake.org/) (version 3.15 or higher)
* [**libsndfile**](http://www.mega-nerd.com/libsndfile/): For loading audio files.
* [**FFTW3**](https://www.fftw.org/): For high-performance Fast Fourier Transforms.

### Build Steps

1.  Clone the repository:
    ```bash
    git clone [https://github.com/Corredor1230/PeriodSplitter.git](https://github.com/Corredor1230/PeriodSplitter.git)
    cd Sitrano
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

Currently, the analysis is configured, compiled, and run from the `main()` function.

1.  **Modify `main.cpp`:**
    * Set the `filename` variable by uncommenting `std::string filename = openFileDialog();` or hard-coding a path.
    * Adjust the parameters in the `AnalysisConfig` and its sub-structs to suit your audio source.

2.  **Compile and Run:**
    ```bash
    cd build
    make
    ./Sitrano
    ```

---

## Output

The program saves the analysis data to three binary files:

* **`INDEX_FILENAME.bin`**: A `std::vector<int>` containing the sample index for the start of each analysis frame. The size of this vector defines the time dimension of the other files.
* **`AMP_FILENAME.bin`**: A `std::vector<std::vector<float>>` (Harmonic x Time) containing the tracked amplitude for each partial in each frame.
* **`FREQ_FILENAME.bin`**: A `std::vector<std::vector<float>>` (Harmonic x Time) containing the tracked frequency for each partial in each frame.

These files can be loaded by other C++ applications or read in scripting languages like Python (using `numpy.fromfile`) for plotting and further modeling.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.