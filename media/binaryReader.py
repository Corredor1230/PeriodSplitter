import argparse
import numpy as np
import matplotlib.pyplot as plt
import soundfile as sf
import os
import pathlib  # Added for easier path handling
import sys      # Added for exit codes

# =============================
# Binary Loader for .sihat
# =============================

def load_sihat_data(path):
    """
    Loads all data from a single consolidated .sihat file.
    
    The file structure is expected to be:
    1. Indices data
    2. Fast RMS data
    3. Slow Frequency data
    """
    with open(path, "rb") as f:
        # --- Block 1: Load Indices ---
        # Read n (number of indices)
        n = np.frombuffer(f.read(4), dtype=np.uint32)[0]
        # Read n * uint32_t
        indices = np.frombuffer(f.read(n * 4), dtype=np.uint32)

        # --- Block 2: Load Fast RMS Data ---
        # Read number of harmonics
        num_harmonics_fast = np.frombuffer(f.read(4), dtype=np.uint32)[0]
        fastData = []
        for _ in range(num_harmonics_fast):
            # Read length of this harmonic's data
            length = np.frombuffer(f.read(4), dtype=np.uint32)[0]
            # Read length * float
            data = np.frombuffer(f.read(length * 4), dtype=np.float32)
            fastData.append(data)

        # --- Block 3: Load Slow Frequency Data ---
        # Read number of harmonics (should be same as fast)
        num_harmonics_slow = np.frombuffer(f.read(4), dtype=np.uint32)[0]
        slowData = []
        
        # Define the C-style struct: { uint32_t index, float value, float ratio }
        # 4 bytes + 4 bytes + 4 bytes = 12 bytes
        cp_dtype = np.dtype([("index", np.uint32), ("value", np.float32), ("ratio", np.float32)])
        
        for _ in range(num_harmonics_slow):
            # Read number of change points for this harmonic
            length = np.frombuffer(f.read(4), dtype=np.uint32)[0]
            # Read length * 12 bytes
            entries = np.frombuffer(f.read(length * cp_dtype.itemsize), dtype=cp_dtype)
            slowData.append(entries)

    print(f"Loaded {len(indices)} indices, {len(fastData)} fast harmonics, and {len(slowData)} slow harmonics.")
    return indices, fastData, slowData

# =============================
# Plotting (Unchanged)
# =============================

def plot_data(indices, fastData, slowData):
    """
    Visualizes:
      - RMS (amplitude) evolution per harmonic
      - Frequency (absolute) evolution per harmonic
      - Frequency ratio (to fundamental) per harmonic
    """

    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)

    # ---------------------------
    # 1️⃣ RMS per harmonic
    # ---------------------------
    for i, rms in enumerate(fastData):
        # Ensure we don't plot more data than we have index labels for
        plot_len = min(len(indices), len(rms))
        axes[0].plot(indices[:plot_len], rms[:plot_len], label=f"H{i+1}")
    axes[0].set_title("Fast RMS (per harmonic)")
    axes[0].set_ylabel("Amplitude (RMS)")
    axes[0].legend(loc="upper right", ncol=4, fontsize=8)
    axes[0].grid(True, alpha=0.3)

    # ---------------------------
    # 2️⃣ Frequency per harmonic
    # ---------------------------
    for i, freqPoints in enumerate(slowData):
        if len(freqPoints) == 0: continue
        freq = np.zeros_like(indices, dtype=float)
        
        for j, point in enumerate(freqPoints):
            idx = point["index"]
            val = point["value"]
            
            # Find the end index for this segment
            is_last_point = (j == len(freqPoints) - 1)
            next_idx = indices[-1] + 1 if is_last_point else freqPoints[j+1]["index"]
            
            mask = (indices >= idx) & (indices < next_idx)
            freq[mask] = val
            
        axes[1].plot(indices, freq, label=f"H{i+1}")
        
    axes[1].set_yscale('log')
    axes[1].set_title("Slow Frequency (Hz)")
    axes[1].set_ylabel("Frequency (Hz)")
    axes[1].legend(loc="upper right", ncol=4, fontsize=8)
    axes[1].grid(True, alpha=0.3)

    # ---------------------------
    # 3️⃣ Harmonic Ratio (f / f₀)
    # ---------------------------
    for i, freqPoints in enumerate(slowData):
        if len(freqPoints) == 0: continue
        ratio = np.zeros_like(indices, dtype=float)
        
        for j, point in enumerate(freqPoints):
            idx = point["index"]
            r = point["ratio"]
            
            is_last_point = (j == len(freqPoints) - 1)
            next_idx = indices[-1] + 1 if is_last_point else freqPoints[j+1]["index"]

            mask = (indices >= idx) & (indices < next_idx)
            ratio[mask] = r
            
        axes[2].plot(indices, ratio, label=f"H{i+1}")
        
    axes[2].set_title("Harmonic Ratio (f / f₀)")
    axes[2].set_xlabel("Sample index (window start)")
    axes[2].set_ylabel("Ratio")
    axes[2].legend(loc="upper right", ncol=4, fontsize=8)
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


# =============================
# Synth mode (Signature Updated)
# =============================
def synthesize(indices, fastData, slowData, sr=48000, filepath="gen/output.wav"):
    """
    Synthesizes the harmonic data into a WAV file.
    
    Args:
        indices (np.array): Array of sample indices.
        fastData (list): List of np.arrays (RMS amplitude envelopes).
        slowData (list): List of structured np.arrays (frequency change points).
        sr (int): The sample rate for the output audio.
        filepath (str): The full path (including directory) to save the .wav file.
    """
    
    if not indices.any():
        print("Error: No indices found. Cannot synthesize.")
        return

    # Create the output directory if it doesn't exist
    saveDir = os.path.dirname(filepath)
    if saveDir:
        os.makedirs(saveDir, exist_ok=True)

    total_samples = indices[-1]
    output = np.zeros(total_samples, dtype=np.float32)

    num_harmonics = min(len(fastData), len(slowData))
    if num_harmonics == 0:
        print("Warning: No harmonic data found to synthesize.")
        return

    for i in range(num_harmonics):
        rms = fastData[i]
        freq_points = slowData[i]
        
        if len(rms) == 0 or len(freq_points) == 0: 
            print(f"Skipping harmonic {i+1}: missing data.")
            continue

        # 1. Create full-length amplitude envelope
        # Interpolate from sparse RMS data to match the indices' length
        amp_at_indices = np.interp(np.arange(len(indices)), np.linspace(0, len(indices)-1, len(rms)), rms)
        # Interpolate from indices to match total sample length
        amp_full = np.interp(np.arange(total_samples), indices, amp_at_indices, left=0.0)

        # 2. Create full-length frequency envelope
        freq_at_indices = np.zeros_like(indices, dtype=float)
        for j, point in enumerate(freq_points):
            idx = point["index"]
            val = point["value"]
            
            is_last_point = (j == len(freq_points)-1)
            next_idx = indices[-1] + 1 if is_last_point else freq_points[j+1]["index"]
            
            # Interpolate frequency between change points
            start_val = val
            end_val = val if is_last_point else freq_points[j+1]["value"]
            
            mask = (indices >= idx) & (indices < next_idx)
            if mask.any():
                interp_freq = np.linspace(start_val, end_val, mask.sum(), endpoint=False)
                freq_at_indices[mask] = interp_freq

        # Interpolate from indices to match total sample length
        freq_full = np.interp(np.arange(total_samples), indices, freq_at_indices)

        # 3. Synthesize oscillator
        # phase = np.cumsum(2 * np.pi * freq_full / sr)
        # More stable phase calculation:
        phase = np.zeros_like(output)
        phase[1:] = np.cumsum(2 * np.pi * freq_full[1:] / sr)
        
        output += amp_full * np.sin(phase)

    # Normalize output
    output /= num_harmonics # Average contribution
    
    target_peak = 0.9
    peak = np.max(np.abs(output))
    if peak > 1e-6: # Avoid division by zero
        output = (output / peak) * target_peak
    else:
        print("Warning: Synthesized audio is silent.")

    sf.write(filepath, output, sr)
    print(f"Successfully synthesized audio to: {filepath}")

# =============================
# Main
# =============================
if __name__ == "__main__":
    inPrefix = "DATA_"
    instr = "AGuit"
    string = "6"
    freq = "82"
    dyn = "F"
    inExt = ".sihat"
    outExt = ".wav"
    outputDir = "gen/"
    inputDir = "sihat/"
    outPrefix = "GEN_"
    sampleRate = 96000

    runMode = "synth"

    fileName = instr + string + "_" + freq + "_" + dyn
    inName = inPrefix + fileName + inExt
    inPath = inputDir + inName
    outName = outPrefix + fileName + outExt
    outPath = outputDir + outName

    if runMode == "synthesize_all":
        outPath = outputDir

    parser = argparse.ArgumentParser(
        description="Read, plot, or synthesize .sihat harmonic data files."
    )
    
    # Renamed 'input_file' to 'input_path'
    parser.add_argument(
        "input_path", 
        help="Path to a .sihat file (for plot/synth modes) or a directory (for synthesize_all mode).",
        # Example default path based on your C++ code
        default=inPath,
        nargs="?" # Makes the default work even if no arg is provided
    )
    
    parser.add_argument(
        "--mode", 
        choices=["plot", "synth", "synthesize_all"], # Added 'synthesize_all'
        default=runMode,
        help="Action to perform. 'synthesize_all' will process all .sihat files in the input_path directory."
    )
    
    parser.add_argument(
        "--sr",
        type=int,
        default=sampleRate,
        help="Sample rate for synthesis."
    )
    
    parser.add_argument(
        "--out",
        type=str,
        default=outPath,
        help="Output file path (for 'synth' mode) or output directory (for 'synthesize_all' mode)."
    )
    
    args = parser.parse_args()
    
    # Handle case where default path might be empty
    input_path_arg = args.input_path if args.input_path else "."
    input_path = pathlib.Path(input_path_arg)


    # --- Main Logic ---
    try:
        if args.mode == "synthesize_all":
            if not input_path.is_dir():
                print(f"Error: For 'synthesize_all' mode, input_path must be a directory.")
                print(f"Provided path: {input_path}")
                sys.exit(1) # Exit with an error code
            
            # Ensure output is a directory
            output_dir = pathlib.Path(args.out)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Find all .sihat files in the directory (non-recursively)
            # Use .rglob("*.sihat") for recursive search (subdirectories)
            sihat_files = list(input_path.glob("*.sihat"))
            
            if not sihat_files:
                print(f"No .sihat files found in directory: {input_path}")
                sys.exit(0) # Not an error, just nothing to do
            
            print(f"Found {len(sihat_files)} files. Starting batch synthesis in '{output_dir}'...")
            
            for file_path in sihat_files:
                print(f"--- Processing: {file_path.name} ---")
                try:
                    indices, fastData, slowData = load_sihat_data(file_path)
                    
                    # Create unique output filename based on input
                    rawName = file_path.stem + ".wav"
                    if rawName.startswith(inPrefix):
                        output_filename = rawName[len(inPrefix):]
                    else:
                        output_filename = rawName
                    output_filepath = output_dir / output_filename
                    
                    synthesize(indices, fastData, slowData, sr=args.sr, filepath=str(output_filepath))
                except Exception as e:
                    # Log error for this file and continue with the next
                    print(f"Failed to process {file_path.name}: {e}")
            print("--- Batch synthesis complete. ---")

        elif args.mode in ["plot", "synth"]:
            if not input_path.is_file():
                print(f"Error: For 'plot' or 'synth' mode, input_path must be a .sihat file.")
                print(f"Provided path: {input_path}")
                sys.exit(1)
                
            print(f"Loading data from: {input_path}")
            indices, fastData, slowData = load_sihat_data(input_path)
            print("Data loaded successfully.")

            if args.mode == "plot":
                plot_data(indices, fastData, slowData)
            else: # mode == "synth"
                # For single synth, args.out is the full filepath
                synthesize(indices, fastData, slowData, sr=args.sr, filepath=args.out)
        
        else:
             print(f"Error: Unknown mode '{args.mode}'")
             sys.exit(1)

    except FileNotFoundError:
        print(f"Error: Input path not found at '{input_path}'")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

