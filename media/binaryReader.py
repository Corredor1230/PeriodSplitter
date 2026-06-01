import struct
import os
from dataclasses import dataclass, field
from typing import List, Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional
import tkinter as tk
from tkinter import filedialog

# ==========================================
# Struct Equivalents (Dataclasses)
# ==========================================

@dataclass
class AnalysisType:
    transient_analysis: bool = False
    harmonic_analysis: bool = False

@dataclass
class Header:
    type: AnalysisType = field(default_factory=AnalysisType)
    sample_rate: int = 0
    f0: float = 0.0
    filename: str = ""

@dataclass
class Harmonic:
    num_frames: int = 0
    indices: List[int] = field(default_factory=list)
    amp: List[List[float]] = field(default_factory=list)
    pha: List[List[float]] = field(default_factory=list)
    freqRate: List[List[float]] = field(default_factory=list)
    rms: float = 0.0

@dataclass
class TransientMeta:
    t_start: int = 0
    t_end: int = 0
    env_hop_size: int = 0
    spec_hop_size: int = 0
    floor_hop_size: int = 0
    spec_window_size: int = 0
    spec_num_bins: int = 0
    spec_num_frames: int = 0
    spec_nfft: int = 0
    num_bands: int = 0
    harm_hop_size: int = 0
    harm_start_sample: int = 0

@dataclass
class Envelope1D:
    first_index: int = 0
    env: List[float] = field(default_factory=list)

@dataclass
class SpectralBin:
    freq: float = 0.0
    mag: float = 0.0
    amp: float = 0.0
    phase: float = 0.0

@dataclass
class EnvelopePoint:
    freqRate: float = 0.0
    crest_factor: float = 0.0
    amp: float = 0.0

@dataclass
class Overtone:
    target: SpectralBin = field(default_factory=SpectralBin)
    envelope: List[EnvelopePoint] = field(default_factory=list)

@dataclass
class Transient:
    meta: TransientMeta = field(default_factory=TransientMeta)
    rise_time: int = 0
    peak_amp: float = 0.0
    envelope: Envelope1D = field(default_factory=Envelope1D)
    centroid: List[float] = field(default_factory=list)
    flatness: List[float] = field(default_factory=list)
    bands: List[float] = field(default_factory=list)
    overtones: List[Overtone] = field(default_factory=list)
    floors: List[float] = field(default_factory=list)
    rms: float = 0.0

# Loader Function
def load_sihat_file(fullpath: str) -> Tuple[Optional[Header], Optional[Harmonic], Optional[Transient]]:
    """
    Reads a binary .sihat file and returns the Header, Harmonic, and Transient data objects.
    Returns (None, None, None) if the file cannot be opened.
    """
    if not os.path.exists(fullpath):
        print(f"Failed to open file for reading: {fullpath}")
        return None, None, None

    header = Header()
    harmonic = Harmonic()
    transient = Transient()

    # Helper function to unpack data safely (using little-endian '<')
    def read_format(f, fmt):
        size = struct.calcsize(fmt)
        data = f.read(size)
        if not data or len(data) != size:
            raise EOFError("Unexpected end of file")
        return struct.unpack(fmt, data)

    # Replaces the C++ lambda readFloatVector
    def read_float_vector(f) -> List[float]:
        size = read_format(f, '<I')[0]
        if size > 0:
            return list(read_format(f, f'<{size}f'))
        return []

    try:
        with open(fullpath, 'rb') as inFile:
            
            # -----------------------------------------
            # HEADER LOAD
            # -----------------------------------------
            t_a, h_a = read_format(inFile, '<BB')
            header.type.transient_analysis = bool(t_a)
            header.type.harmonic_analysis = bool(h_a)

            header.sample_rate, header.f0, string_size = read_format(inFile, '<IfI')
            
            if string_size > 0:
                header.filename = inFile.read(string_size).decode('utf-8', errors='ignore')

            # -----------------------------------------
            # HARMONIC ANALYSIS
            # -----------------------------------------
            if header.type.harmonic_analysis:
                # Block 1: Indices
                harmonic.num_frames = read_format(inFile, '<I')[0]
                if harmonic.num_frames > 0:
                    harmonic.indices = list(read_format(inFile, f'<{harmonic.num_frames}I'))

                # Block 2, 3, 4 helper for 2D vectors (Amp, Pha, Freq)
                def read_2d_float_vector(f):
                    num_harmonics = read_format(f, '<I')[0]
                    result = []
                    for _ in range(num_harmonics):
                        num_frames = read_format(f, '<I')[0]
                        if num_frames > 0:
                            result.append(list(read_format(f, f'<{num_frames}f')))
                        else:
                            result.append([])
                    return result

                # Block 2: Amp
                harmonic.amp = read_2d_float_vector(inFile)
                # Block 3: Phase
                harmonic.pha = read_2d_float_vector(inFile)
                # Block 4: Frequency
                harmonic.freqRate = read_2d_float_vector(inFile)

            # -----------------------------------------
            # TRANSIENT ANALYSIS
            # -----------------------------------------
            if header.type.transient_analysis:
                # 1. Read the Range and Metadata
                (
                    transient.meta.t_start, transient.meta.t_end,
                    transient.meta.env_hop_size, transient.meta.spec_hop_size,
                    transient.meta.floor_hop_size, transient.meta.spec_window_size,
                    transient.meta.spec_num_bins, transient.meta.spec_num_frames,
                    transient.meta.spec_nfft
                ) = read_format(inFile, '<IIIIIIIII')

                # 2. Read scalar metrics
                transient.rise_time, transient.peak_amp = read_format(inFile, '<if')

                # 3. Read 1D analysis envelopes
                env_size = read_format(inFile, '<I')[0]
                if env_size > 0:
                    transient.envelope.first_index = read_format(inFile, '<I')[0]
                    transient.envelope.env = list(read_format(inFile, f'<{env_size}f'))

                transient.centroid = read_float_vector(inFile)
                transient.flatness = read_float_vector(inFile)

                # 4. Read the Band Partials
                transient.meta.num_bands = read_format(inFile, '<I')[0]

                # 5. Read the flattened band envelopes
                transient.bands = read_float_vector(inFile)

                # Metadata for the harmonic section
                transient.meta.harm_hop_size, transient.meta.harm_start_sample = read_format(inFile, '<II')

                num_overtones = read_format(inFile, '<I')[0]
                for _ in range(num_overtones):
                    overtone = Overtone()
                    # target (SpectralBin)
                    (
                        overtone.target.freq, overtone.target.mag, 
                        overtone.target.amp, overtone.target.phase
                    ) = read_format(inFile, '<dddd')

                    # envelope (EnvelopePoint)
                    frame_num = read_format(inFile, '<I')[0]
                    for _ in range(frame_num):
                        f_val, c_val, a_val = read_format(inFile, '<fff')
                        overtone.envelope.append(EnvelopePoint(freqRate=f_val, crest_factor=c_val, amp=a_val))
                    
                    transient.overtones.append(overtone)

                # Floor values
                transient.floors = read_float_vector(inFile)

            # -----------------------------------------
            # RMS BLOCK
            # -----------------------------------------
            if header.type.harmonic_analysis:
                harmonic.rms = read_format(inFile, '<f')[0]
            
            if header.type.transient_analysis:
                transient.rms = read_format(inFile, '<f')[0]

    except EOFError as e:
        print(f"File parsing error (EOF): {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    return header, harmonic, transient

def plot_sihat_feature(header, harmonic, transient, feature: str):
    """
    Plots specific features from the Sihat data structures.
    Creates a new figure for each call so multiple plots can be generated consecutively.
    """
    # Safety check
    if header is None:
        print("Invalid data: Header is None. Cannot plot.")
        return

    plt.figure(figsize=(10, 6))
    
    # ---------------------------------------------------------
    # 1. Harmonic Frequency (har_freq)
    # ---------------------------------------------------------
    if feature == 'har_freq':
        if not header.type.harmonic_analysis or not harmonic:
            print("No harmonic data available for 'har_freq'.")
            return
            
        plt.title(f"Harmonic Frequencies - {header.filename}")
        plt.xlabel("Time (Frame Index or Sample Index)")
        plt.ylabel("Frequency (Hz)")
        
        # Filter settings
        db_threshold = -60.0
        consecutive_frames_limit = 10
        
        for i, freqs in enumerate(harmonic.freq):
            amps = harmonic.amp[i]
            
            # Default cutoff is the full length of the overtone
            cutoff_idx = len(freqs)
            
            # Safety check: ensure we have matching amplitude data
            if len(amps) == len(freqs):
                # Convert to numpy array and calculate dB (avoiding log10(0))
                amps_array = np.array(amps)
                amps_db = 20 * np.log10(np.maximum(amps_array, 1e-10))
                
                below_threshold_count = 0
                
                for j, db_val in enumerate(amps_db):
                    if db_val < db_threshold:
                        below_threshold_count += 1
                        if below_threshold_count >= consecutive_frames_limit:
                            # Set cutoff to exactly where the 10-frame drop started
                            cutoff_idx = j - consecutive_frames_limit + 1
                            break
                    else:
                        below_threshold_count = 0
            else:
                print(f"Warning: Amp and Freq length mismatch for overtone {i}. Skipping filter.")

            # Slice the data up to the calculated cutoff index
            plot_freqs = freqs[:cutoff_idx]
            
            # If lengths match the global indices array, use it for X.
            if len(freqs) == len(harmonic.indices):
                plot_indices = harmonic.indices[:cutoff_idx]
                plt.plot(plot_indices, plot_freqs, label=f'Overtone {i}')
            else:
                plt.plot(plot_freqs, label=f'Overtone {i}')
                
    # ---------------------------------------------------------
    # 2. Harmonic Amplitude (har_amp)
    # ---------------------------------------------------------
    elif feature == 'har_amp':
        if not header.type.harmonic_analysis or not harmonic:
            print("No harmonic data available for 'har_amp'.")
            return
            
        plt.title(f"Harmonic Amplitudes - {header.filename}")
        plt.xlabel("Time (Frame Index or Sample Index)")
        plt.ylabel("Amplitude")
        
        for i, amps in enumerate(harmonic.amp):
            if len(amps) == len(harmonic.indices):
                plt.plot(harmonic.indices, amps, label=f'Overtone {i}')
            else:
                plt.plot(amps, label=f'Overtone {i}')

    # ---------------------------------------------------------
    # 3. Amplitude Envelope (amp_env)
    # ---------------------------------------------------------
    elif feature == 'amp_env':
        if not header.type.transient_analysis or not transient:
            print("No transient data available for 'amp_env'.")
            return
            
        plt.title(f"Transient Amplitude Envelope - {header.filename}")
        plt.xlabel("Time (Samples)")
        plt.ylabel("Amplitude")
        
        env = transient.envelope.env
        t_start = transient.meta.t_start
        hop_size = transient.meta.env_hop_size or 1 # Fallback to 1 to avoid div zero
        
        # X-axis array starting at tStart
        x_axis = np.arange(len(env)) * hop_size + t_start
        plt.plot(x_axis, env, label='Envelope')
        
        # Red vertical line for Rise Time (measured from tStart)
        peak_time = t_start + transient.rise_time
        plt.axvline(x=peak_time, color='r', linestyle='--', label=f'Peak/Rise Time ({peak_time})')
        plt.legend()

    # ---------------------------------------------------------
    # 4. Spectral Flatness (flatness)
    # ---------------------------------------------------------
    elif feature == 'flatness':
        if not header.type.transient_analysis or not transient:
            print("No transient data available for 'flatness'.")
            return
            
        plt.title(f"Spectral Flatness - {header.filename}")
        plt.xlabel("Time (Samples)")
        plt.ylabel("Flatness")
        
        flatness = transient.flatness
        t_start = transient.meta.t_start
        hop_size = transient.meta.spec_hop_size or 1
        
        x_axis = np.arange(len(flatness)) * hop_size + t_start
        plt.plot(x_axis, flatness, color='orange')

    # ---------------------------------------------------------
    # 5. Spectral Centroid (centroid)
    # ---------------------------------------------------------
    elif feature == 'centroid':
        if not header.type.transient_analysis or not transient:
            print("No transient data available for 'centroid'.")
            return
            
        plt.title(f"Spectral Centroid - {header.filename}")
        plt.xlabel("Time (Samples)")
        plt.ylabel("Centroid (Hz)")
        
        centroid = transient.centroid
        t_start = transient.meta.t_start
        hop_size = transient.meta.spec_hop_size or 1
        
        x_axis = np.arange(len(centroid)) * hop_size + t_start
        plt.plot(x_axis, centroid, color='green')

    # ---------------------------------------------------------
    # 6. Band Envelopes (band_env)
    # ---------------------------------------------------------
    elif feature == 'band_env':
        if not header.type.transient_analysis or not transient:
            print("No transient data available for 'band_env'.")
            return
            
        plt.title(f"Band Envelopes - {header.filename}")
        plt.xlabel("Time (Samples)")
        plt.ylabel("Amplitude")
        
        bands = transient.bands
        num_bands = transient.meta.num_bands
        t_start = transient.meta.t_start
        hop_size = transient.meta.env_hop_size or 1
        
        if num_bands > 0 and len(bands) > 0:
            frames_per_band = len(bands) // num_bands
            x_axis = np.arange(frames_per_band) * hop_size + t_start
            
            # Unflatten and plot
            for i in range(num_bands):
                start_idx = i * frames_per_band
                end_idx = start_idx + frames_per_band
                plt.plot(x_axis, bands[start_idx:end_idx], label=f'Band {i}')
        else:
            print("No band data to plot.")

    # ---------------------------------------------------------
    # 7. Transient Harmonic Frequencies (t_har_freq)
    # ---------------------------------------------------------
    elif feature == 't_har_freq':
        if not header.type.transient_analysis or not transient:
            print("No transient data available for 't_har_freq'.")
            return
            
        plt.title(f"Transient Harmonic Frequencies - {header.filename}")
        plt.xlabel("Time (Samples)")
        plt.ylabel("Frequency (Hz)")
        
        t_start = transient.meta.harm_start_sample
        hop_size = transient.meta.harm_hop_size or 1
        
        for i, overtone in enumerate(transient.overtones):
            env_freqs = [pt.freq for pt in overtone.envelope]
            x_axis = np.arange(len(env_freqs)) * hop_size + t_start
            plt.plot(x_axis, env_freqs, label=f'Transient Overtone {i}')
    
    # ---------------------------------------------------------
    # 8. Transient Harmonic Amplitudes (t_har_amp)
    # ---------------------------------------------------------
    elif feature == 't_har_amp':
        if not header.type.transient_analysis or not transient:
            print("No transient data available for 't_har_amp'.")
            return
            
        plt.title(f"Transient Harmonic Amplitudes - {header.filename}")
        plt.xlabel("Time (Samples)")
        plt.ylabel("Amplitude")
        
        t_start = transient.meta.harm_start_sample
        hop_size = transient.meta.harm_hop_size or 1
        
        for i, overtone in enumerate(transient.overtones):
            env_amps = [pt.amp for pt in overtone.envelope]
            x_axis = np.arange(len(env_amps)) * hop_size + t_start
            plt.plot(x_axis, env_amps, label=f'Transient Overtone {i}')
            
    else:
        print(f"Unknown feature requested: {feature}")
        plt.close() # Close the figure if the feature is invalid

    plt.grid(True, alpha=0.3)
    plt.tight_layout()

def get_filename():
    root = tk.Tk()
    root.withdraw()
    
    # 2. Open the file explorer dialog
    print("Waiting for file selection...")
    filepath = filedialog.askopenfilename(
        title="Select a Sihat Audio File",
        filetypes=[("Sihat Files", "*.sihat"), ("All Files", "*.*")]
    )
    
    # If the user clicks "Cancel" or closes the window, filepath will be empty
    if not filepath:
        print("No file selected. Exiting.")
        return

    return filepath

# Main Execution Example
def main():
    # 1. Define the path to your binary file
    filepath = get_filename()
    
    # 2. Load the data using the previously written function
    header, harmonic, transient = load_sihat_file(filepath)
    if header is None:
        print("Failed to load file, exiting.")
        return

    # 3. Request your plots consecutively
    print("Generating plots...")
    
    # Allowed: har_freq, har_amp, amp_env, flatness, centroid, band_env, t_har_freq
    plot_sihat_feature(header, harmonic, transient, 'har_freq')
    plot_sihat_feature(header, harmonic, transient, 'har_amp')
    plot_sihat_feature(header, harmonic, transient, 'amp_env')
    plot_sihat_feature(header, harmonic, transient, 'band_env')
    
    # 4. Display all requested plots at once
    plt.show()

if __name__ == "__main__":
    main()