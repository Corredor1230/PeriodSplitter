import argparse
import numpy as np
import matplotlib.pyplot as plt
import soundfile as sf
import os
import pathlib
import sys
from scipy.interpolate import interp1d

# Binary Loader for .sihat
def load_sihat_data(path):
    """
    Loads all data from a single consolidated .sihat file.
    Structure:
    0. Fundamental frequency (f0)
    1. Indices data
    2. Harmonic Amplitude data (Fast)
    3. Harmonic Frequency data (Slow)
    4. Transient Data (Range + Scalogram)
    """
    with open(path, "rb") as f:
        # Block 0: Load f0
        f0 = np.frombuffer(f.read(4), dtype=np.float32)[0]

        # Block 1: Load Indices
        n = np.frombuffer(f.read(4), dtype=np.uint32)[0]
        indices = np.frombuffer(f.read(n * 4), dtype=np.uint32)

        # Block 2: Load RMS Data
        num_harmonics_fast = np.frombuffer(f.read(4), dtype=np.uint32)[0]
        ampData = []
        for _ in range(num_harmonics_fast):
            length = np.frombuffer(f.read(4), dtype=np.uint32)[0]
            data = np.frombuffer(f.read(length * 4), dtype=np.float32)
            ampData.append(data)

        # Block 3: Load Frequency Data
        num_harmonics_slow = np.frombuffer(f.read(4), dtype=np.uint32)[0]
        freqData = []
        ratios = []
        cp_dtype = np.dtype([("index", np.uint32), ("value", np.float32), ("ratio", np.float32)])
        
        for _ in range(num_harmonics_slow):
            length = np.frombuffer(f.read(4), dtype=np.uint32)[0]
            if length > 0:
                entries = np.frombuffer(f.read(length * cp_dtype.itemsize), dtype=cp_dtype)
                freqData.append(entries)
                ratios.append(entries['ratio'])
            else:
                freqData.append(np.array([], dtype=cp_dtype))

        # Block 4: Load Transient Data
        # Read Range
        t_start = np.frombuffer(f.read(4), dtype=np.uint32)[0]
        t_end = np.frombuffer(f.read(4), dtype=np.uint32)[0]
        
        # Read Scalogram Partials
        num_partials = np.frombuffer(f.read(4), dtype=np.uint32)[0]
        transient_partials = []

        for _ in range(num_partials):
            freq = np.frombuffer(f.read(4), dtype=np.float32)[0]
            hop = np.frombuffer(f.read(4), dtype=np.uint32)[0]
            d_len = np.frombuffer(f.read(4), dtype=np.uint32)[0]
            
            data = np.array([], dtype=np.float32)
            if d_len > 0:
                data = np.frombuffer(f.read(d_len * 4), dtype=np.float32)
            
            transient_partials.append({
                "freq": freq,
                "hop": hop,
                "data": data
            })
        
        hRms = np.frombuffer(f.read(4), dtype=np.float32)[0]
        tRms = np.frombuffer(f.read(4), dtype=np.float32)[0]

        rmsData = {
            "hRms": hRms,
            "tRms": tRms
        }

        transientData = {
            "start": t_start,
            "end": t_end,
            "partials": transient_partials
        }

    print(f"Loaded {len(indices)} indices, {len(ampData)} fast harmonics.")
    print(f"Transient found: {t_start} to {t_end} with {len(transient_partials)} partials.")
    return indices, ampData, freqData, f0, ratios, transientData, rmsData

# Plotting
def plot_data(indices, fastData, slowData, transientData, 
              plots=["rms", "ratio", "scalogram"], 
              stacked=True,
              amp_threshold=0.0005,
              sr=96000): # <--- Added SR for Time conversion
    
    if not plots:
        print("No plots requested.")
        return

    # ----------------------------------------------------
    # 1. ACADEMIC STYLE SETTINGS
    # ----------------------------------------------------
    plt.rcParams.update({
        'font.size': 24,
        'axes.titlesize': 24,
        'axes.labelsize': 24,
        'xtick.labelsize': 24,
        'ytick.labelsize': 24,
        'legend.fontsize': 24,
        'lines.linewidth': 2,
        'figure.titlesize': 30,
        'font.family': 'sans-serif'
    })

    # Helper to convert samples to seconds
    def to_sec(idx_array):
        return np.array(idx_array) / sr

    # ----------------------------------------------------
    # 2. SETUP PLOTS
    # ----------------------------------------------------
    axes_list = []
    if stacked:
        num_plots = len(plots)
        fig, axes = plt.subplots(num_plots, 1, figsize=(12, 4 * num_plots), sharex=False)
        if num_plots == 1: axes_list = [axes]
        else: axes_list = list(axes)
        # Optional: Main title often removed for papers
        # fig.suptitle(f"Analysis Data (Threshold: {amp_threshold})")

    for i, plot_type in enumerate(plots):
        if stacked:
            ax = axes_list[i]
        else:
            fig = plt.figure(figsize=(10, 5))
            ax = plt.gca()
            fig.canvas.manager.set_window_title(f"Detailed View: {plot_type.upper()}")

        # -----------------------------
        # RMS PLOT
        # -----------------------------
        if plot_type == "rms":
            # Convert X-axis (indices) to Seconds
            time_axis = to_sec(indices)
            
            for k, rms in enumerate(fastData):
                plot_len = min(len(time_axis), len(rms))
                ax.plot(time_axis[:plot_len], rms[:plot_len], label=f"H{k+1}")
            
            ax.set_title("Spectral Envelopes", fontweight='bold')
            ax.set_ylabel("Amplitude")
            ax.set_xlabel("Time (s)") # Explicit label
            ax.grid(True, alpha=0.3)

        # -----------------------------
        # FREQUENCY & RATIO PLOTS
        # -----------------------------
        elif plot_type in ["freq", "ratio"]:
            data_key = "value" if plot_type == "freq" else "ratio"
            
            # Convert X-axis (indices) to Seconds
            time_axis = to_sec(indices)

            for k, freqPoints in enumerate(slowData):
                if len(freqPoints) == 0: continue
                
                # Reconstruct Y values
                y_values = np.zeros_like(indices, dtype=float)
                
                for j, point in enumerate(freqPoints):
                    idx = point["index"]
                    val = point[data_key]
                    
                    is_last_point = (j == len(freqPoints) - 1)
                    next_idx = indices[-1] + 1 if is_last_point else freqPoints[j+1]["index"]
                    
                    mask = (indices >= idx) & (indices < next_idx)
                    y_values[mask] = val
                
                # Noise Gate
                if k < len(fastData):
                    rms = fastData[k]
                    limit = min(len(y_values), len(rms))
                    silence_mask = rms[:limit] < amp_threshold
                    y_values[:limit][silence_mask] = np.nan

                ax.plot(time_axis, y_values, label=f"H{k+1}")

            if plot_type == "freq":
                ax.set_yscale('log')
                ax.set_title("Harmonic Frequencies", fontweight='bold')
                ax.set_ylabel("Frequency (Hz)")
            else:
                ax.set_title("Harmonic Ratio (f / f0)", fontweight='bold')
                ax.set_ylabel("Ratio")
            
            ax.set_xlabel("Time (s)")
            ax.grid(True, alpha=0.3)

        # -----------------------------
        # SCALOGRAM (Amplitude View)
        # -----------------------------
        elif plot_type == "scalogram":
            t_start = transientData["start"]
            t_end = transientData["end"]
            
            for p in transientData["partials"]:
                if len(p["data"]) == 0: continue
                
                # Create time vector in Seconds
                # Start Time (s) -> End Time (s)
                t_sec_start = t_start / sr
                t_sec_end = t_end / sr
                
                time_vec = np.linspace(t_sec_start, t_sec_end, num=len(p["data"]))
                
                ax.plot(time_vec, p["data"], alpha=0.7, linewidth=1.5)

            ax.set_title("Transient Partials (Envelope)", fontweight='bold')
            ax.set_ylabel("Amplitude")
            ax.set_xlabel("Time (s)")
            ax.grid(True, alpha=0.3)

        if not stacked:
            #plt.tight_layout()
            plt.show()

    if stacked:
        #plt.tight_layout()
        plt.show()

# =============================
# Synthesis Logic
# =============================

def synthesize_transient(transientData, sr, total_samples):
    """
    Generates the audio buffer for the transient section.
    """
    output = np.zeros(total_samples, dtype=np.float32)
    
    start = int(transientData["start"])
    end = int(transientData["end"])
    length = end - start
    
    if length <= 0: return output

    # Buffer for just the transient part
    t_buffer = np.zeros(length, dtype=np.float32)
    time_indices = np.arange(length) # 0 to length-1

    for p in transientData["partials"]:
        freq = p["freq"]
        hop = p["hop"]
        env_data = p["data"]
        
        if len(env_data) == 0: continue

        # Upsample envelope: The envelope data is decimated by 'hop'.
        # We need to stretch it to fit the transient length.
        
        # Original indices of the envelope points
        env_x = np.arange(len(env_data)) * hop
        
        # Target indices (every sample in the transient)
        target_x = np.arange(length)
        
        # Interpolate
        # Use fill_value="extrapolate" to handle edges if hop alignment isn't perfect
        interpolator = interp1d(env_x, env_data, kind='linear', bounds_error=False, fill_value=0.0)
        envelope = interpolator(target_x)
        
        # Generate Sine
        # Phase is simple here because freq is constant for the duration of the transient slice
        phase = 2 * np.pi * freq * time_indices / sr
        sine_wave = np.sin(phase)
        
        t_buffer += sine_wave * envelope

    # Inject into main buffer
    # Ensure we don't go out of bounds
    write_end = min(end, total_samples)
    write_len = write_end - start
    
    if write_len > 0:
        output[start:start+write_len] += t_buffer[:write_len]
        
    return output

def synthesize(f0, ratios, indices, fastData, slowData, transientData, 
               unifiedRMS=None,               # <--- Optional RMS struct
               use_transient=True,            # <--- FLAG 1
               use_harmonics=True,            # <--- FLAG 2
               sr=48000, filepath="gen/output.wav"):
    
    if not indices.any() and transientData["end"] == 0:
        print("Error: No data found. Cannot synthesize.")
        return

    # Create directory
    saveDir = os.path.dirname(filepath)
    if saveDir:
        os.makedirs(saveDir, exist_ok=True)

    # 1. Determine Total Length
    max_harmonic_idx = indices[-1] if len(indices) > 0 else 0
    max_transient_idx = transientData["end"]
    total_samples = max(max_harmonic_idx, max_transient_idx) + 4096 
    
    # Initialize main buffer
    final_output = np.zeros(total_samples, dtype=np.float32)

    # ---------------------------------------------------------
    # 2. Synthesize Transient (If Requested)
    # ---------------------------------------------------------
    if use_transient:
        print("Synthesizing Transient...")
        raw_transient = synthesize_transient(transientData, sr, total_samples)

        # Apply RMS Normalization (if UnifiedRMS data exists)
        if unifiedRMS:
            tStart = int(transientData["start"])
            tEnd = int(transientData["end"])
            t_slice = raw_transient[tStart:tEnd]
            gen_t_rms = np.sqrt(np.mean(t_slice**2)) if len(t_slice) > 0 else 0.0
            
            if gen_t_rms > 1e-9 and unifiedRMS["tRms"] > 0:
                gain = unifiedRMS["tRms"] / gen_t_rms
                raw_transient *= gain
                print(f"   -> Transient Gain Applied: {gain:.2f}")

        final_output += raw_transient
    else:
        print("Skipping Transient...")

    # ---------------------------------------------------------
    # 3. Synthesize Harmonics (If Requested)
    # ---------------------------------------------------------
    if use_harmonics:
        print("Synthesizing Harmonics...")
        harmonic_output = np.zeros(total_samples, dtype=np.float32)
        num_harmonics = min(len(fastData), len(slowData))

        if num_harmonics > 0 and len(indices) > 0:
            for i in range(num_harmonics):
                rms = fastData[i]
                freq_points = slowData[i]
                
                if len(rms) == 0 or len(freq_points) == 0: continue

                # Interpolate Amp
                amp_at_indices = np.interp(np.arange(len(indices)), np.linspace(0, len(indices)-1, len(rms)), rms)
                amp_full = np.zeros(total_samples, dtype=np.float32)
                
                start_idx = indices[0]
                end_idx = indices[-1]
                range_len = end_idx - start_idx + 1
                
                if range_len <= 0: continue

                x_src = indices - start_idx
                x_dst = np.arange(range_len)
                amp_segment = np.interp(x_dst, x_src, amp_at_indices)
                amp_full[start_idx : start_idx+range_len] = amp_segment

                # Interpolate Freq
                freq_segment = np.zeros(range_len, dtype=np.float32)
                curr_freq_val = freq_points[0]["ratio"] * f0
                abs_indices_segment = np.arange(start_idx, end_idx + 1)
                freq_segment[:] = curr_freq_val 
                
                for k in range(len(freq_points) - 1):
                    p_curr = freq_points[k]
                    p_next = freq_points[k+1]
                    
                    idx_a = p_curr["index"]
                    idx_b = p_next["index"]
                    val_a = p_curr["ratio"] * f0
                    val_b = p_next["ratio"] * f0
                    
                    rel_a = int(idx_a - start_idx)
                    rel_b = int(idx_b - start_idx)
                    rel_a = max(0, rel_a)
                    rel_b = min(range_len, rel_b)
                    
                    if rel_b > rel_a:
                        ramp = np.linspace(val_a, val_b, rel_b - rel_a, endpoint=False)
                        freq_segment[rel_a : rel_b] = ramp
                
                last_p = freq_points[-1]
                rel_last = int(last_p["index"] - start_idx)
                if rel_last < range_len:
                    freq_segment[rel_last:] = last_p["ratio"] * f0

                # Generate Phase & Sine
                phase_segment = np.cumsum(2 * np.pi * freq_segment / sr)
                sine_segment = np.sin(phase_segment)
                harmonic_output[start_idx : start_idx+range_len] += amp_full[start_idx : start_idx+range_len] * sine_segment

        # Soft Normalization of additive sum
        if num_harmonics > 0:
            harmonic_output /= np.sqrt(num_harmonics) 

        # Apply RMS Normalization (if UnifiedRMS data exists)
        if unifiedRMS:
            hStart = int(transientData["end"])
            hEnd = min(len(harmonic_output), hStart + sr)
            #Slice is of size SR, meaning 1 second
            h_slice = harmonic_output[hStart:hEnd]
            gen_h_rms = np.sqrt(np.mean(h_slice**2)) if len(h_slice) > 0 else 0.0

            if gen_h_rms > 1e-9 and unifiedRMS["hRms"] > 0:
                gain = unifiedRMS["hRms"] / gen_h_rms
                harmonic_output *= gain
                print(f"   -> Harmonic Gain Applied: {gain:.2f}")

        # Mix into final output
        final_output += harmonic_output
    else:
        print("Skipping Harmonics...")

    # ---------------------------------------------------------
    # 4. Final Output Processing
    # ---------------------------------------------------------
    # Global Limiter
    peak = np.max(np.abs(final_output))
    if peak > 1.0:
        final_output = (final_output / peak) * 0.99

    current_peak = np.max(np.abs(final_output))
    
    # Target peak (e.g., -0.5 dB)
    target_peak = 0.95

    if current_peak > 1e-9:
        # Normalize: Boost (or cut) the signal so the loudest sample hits target_peak
        final_output = (final_output / current_peak) * target_peak
        print(f"Global Normalization: Peak {current_peak:.4f} -> {target_peak}")
    else:
        print("Warning: Output is silent.")

    sf.write(filepath, final_output, sr)
    print(f"Successfully synthesized audio to: {filepath}")

# =============================
# Main
# =============================
if __name__ == "__main__":
    inPrefix = "DATA_"
    instr = "EGuit"
    string = "6"
    freq = "82"
    newF0 = 82.0
    dyn = "F"
    inExt = ".sihat"
    outExt = ".wav"
    outputDir = "tests/withTransient/"
    inputDir = "sihat/withTransient/"
    outPrefix = "GEN" + "_"

    plotsStrings=["rms", "ratio", "scalogram"] # Default plots
    stackplots = False
    ampThreshold = 0.001

    sampleRate = 96000
    replaceF0 = False
    if replaceF0:
        outPrefix = "GENto" + str(int(newF0)) + "_"
    useTransient = True
    useOvertones = True

    runMode = "synth"

    fileName = instr + string + "_" + freq + "_" + dyn
    inName = inPrefix + fileName + inExt
    inPath = inputDir + inName
    outName = outPrefix + fileName + outExt
    outPath = outputDir + outName

    if runMode == "synthesize_all":
        inPath = inputDir
        outPath = outputDir

    parser = argparse.ArgumentParser()
    parser.add_argument("input_path", default=inPath, nargs="?")
    parser.add_argument("--mode", choices=["plot", "synth", "synthesize_all"], default=runMode)
    parser.add_argument("--sr", type=int, default=sampleRate)
    parser.add_argument("--out", type=str, default=outPath)
    
    args = parser.parse_args()
    input_path = pathlib.Path(args.input_path if args.input_path else ".")

    try:
        if args.mode == "synthesize_all":
            if not input_path.is_dir():
                print(f"Error: input_path must be a directory for batch mode.")
                sys.exit(1)
            
            output_dir = pathlib.Path(args.out)
            output_dir.mkdir(parents=True, exist_ok=True)
            sihat_files = list(input_path.glob("*.sihat"))
            
            print(f"Found {len(sihat_files)} files. Processing...")
            
            for file_path in sihat_files:
                print(f"--- {file_path.name} ---")
                try:
                    indices, fast, slow, f0, ratios, tData, rmsData = load_sihat_data(file_path)
                    
                    rawName = file_path.stem + ".wav"
                    if rawName.startswith(inPrefix):
                        output_filename = rawName[len(inPrefix):]
                    else:
                        output_filename = rawName
                    
                    target_f0 = newF0 if replaceF0 else f0
                    synthesize(target_f0, ratios, indices, fast, slow, tData, use_harmonics=useOvertones, use_transient=useTransient, sr=args.sr, filepath=str(output_dir / output_filename), unifiedRMS=rmsData)
                except Exception as e:
                    print(f"Failed: {e}")

        elif args.mode in ["plot", "synth"]:
            if not input_path.is_file():
                print(f"Error: input_path must be a file.")
                sys.exit(1)
                
            print(f"Reading: {input_path}")
            indices, fast, slow, f0, ratios, tData, rmsData = load_sihat_data(input_path)
            
            if args.mode == "plot":
                plot_data(indices, fast, slow, tData, plotsStrings, stackplots, ampThreshold)
            else:
                target_f0 = newF0 if replaceF0 else f0
                synthesize(target_f0, ratios, indices, fast, slow, tData, use_harmonics=useOvertones, use_transient=useTransient,  sr=args.sr, filepath=str(args.out), unifiedRMS=rmsData)
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)