import librosa
import librosa.display
import matplotlib.pyplot as plt
import numpy as np
import os

# ==========================================
# CONFIGURATION
# ==========================================

# 1. Comparison Settings
NUM_FILES = 2              # How many files to plot side-by-side
ANALYSIS_MODE = "linear"      # Options: "mel" or "linear"

# 2. Files to Analyze
FILE_PATHS = [
    "sihat/newTests/GENto329_EGuit6_82_F.wav",
    "sihat/newTests/GEN_EGuit1_329_F.wav",
    # ... add more files here
]

# 3. DSP Settings
N_FFT = 2048
HOP_LENGTH = 512
N_MELS = 128      # Only used if mode is "mel"

# ==========================================
# MAIN LOGIC
# ==========================================

def compare_spectrograms(files, limit, mode):
    # 1. Select the files
    selected_files = files[:limit]
    
    if not selected_files:
        print("No files provided or file list is empty.")
        return

    num_plots = len(selected_files)
    
    # 2. Setup Figure
    # We dynamically create subplots based on the limit
    fig, axes = plt.subplots(1, num_plots, figsize=(6 * num_plots, 5), sharey=True)
    
    if num_plots == 1:
        axes = [axes]

    fig.suptitle(f"Comparison: {mode.upper()} Spectrogram", fontsize=16)

    # 3. Iterate and Plot
    for i, file_path in enumerate(selected_files):
        ax = axes[i]
        
        try:
            print(f"Analyzing {i+1}/{num_plots}: {file_path}")
            y, sr = librosa.load(file_path)
            
            # --- MEL SPECTROGRAM (Human Hearing Scale) ---
            if mode == "mel":
                # Compute power spectrogram (amplitude squared)
                S = librosa.feature.melspectrogram(y=y, sr=sr, n_fft=N_FFT, 
                                                 hop_length=HOP_LENGTH, n_mels=N_MELS)
                
                # Convert to dB (log scale)
                S_dB = librosa.power_to_db(S, ref=np.max)
                
                # Plot
                img = librosa.display.specshow(S_dB, x_axis='time', y_axis='mel', sr=sr, 
                                             hop_length=HOP_LENGTH, ax=ax, cmap='magma')
                ax.set_title(f"{os.path.basename(file_path)}\n(Mel Scale)")

            # --- LINEAR SPECTROGRAM (Standard FFT) ---
            elif mode == "linear":
                # Compute Short-Time Fourier Transform (STFT)
                D = librosa.stft(y, n_fft=N_FFT, hop_length=HOP_LENGTH)
                
                # Convert complex values to Magnitude, then to dB
                S_dB = librosa.amplitude_to_db(np.abs(D), ref=np.max)
                
                # Plot
                # y_axis='linear' plots the raw Hz frequencies evenly
                img = librosa.display.specshow(S_dB, x_axis='time', y_axis='linear', sr=sr, 
                                             hop_length=HOP_LENGTH, ax=ax, cmap='magma')
                ax.set_title(f"{os.path.basename(file_path)}\n(Linear FFT)")

            # Add Colorbar (Legend) only on the last plot to keep it clean
            if i == num_plots - 1:
                fig.colorbar(img, ax=ax, format='%+2.0f dB')
            
            ax.set_xlabel("Time")
            if i == 0:
                ax.set_ylabel("Frequency (Hz)")

        except Exception as e:
            ax.text(0.5, 0.5, "File Error", ha='center', va='center')
            print(f"Error processing {file_path}: {e}")

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    compare_spectrograms(FILE_PATHS, NUM_FILES, ANALYSIS_MODE)