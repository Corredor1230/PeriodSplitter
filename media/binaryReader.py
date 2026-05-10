import argparse
import numpy as np
import matplotlib.pyplot as plt
import soundfile as sf
import os
import pathlib
import sys
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt

# =============================
# Helper: Read float vectors
# =============================
def read_float_vector(f):
    size = np.frombuffer(f.read(4), dtype=np.uint32)[0]
    if size > 0:
        return np.frombuffer(f.read(size * 4), dtype=np.float32)
    return np.array([], dtype=np.float32)

# =============================
# Binary Loader for .sihat
# =============================
def load_sihat_data(path, source_separation=False):
    """
    Loads all data from a single consolidated .sihat file.
    """
    with open(path, "rb") as f:
        # Block 0: Load f0
        f0 = np.frombuffer(f.read(4), dtype=np.float32)[0]

        # Block 1: Load Indices
        n = np.frombuffer(f.read(4), dtype=np.uint32)[0]
        indices = np.frombuffer(f.read(n * 4), dtype=np.uint32)

        # Block 2: Load Fast Harmonic Amp Data
        num_harmonics_fast = np.frombuffer(f.read(4), dtype=np.uint32)[0]
        ampData = []
        for _ in range(num_harmonics_fast):
            length = np.frombuffer(f.read(4), dtype=np.uint32)[0]
            data = np.frombuffer(f.read(length * 4), dtype=np.float32)
            ampData.append(data)

        # Block 3: Load Slow Harmonic Freq Data
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
        t_start = np.frombuffer(f.read(4), dtype=np.uint32)[0]
        t_end = np.frombuffer(f.read(4), dtype=np.uint32)[0]

        if source_separation:
            # NEW FORMAT: STransientResults

            env_hop = np.frombuffer(f.read(4), dtype=np.uint32)[0]
            spec_hop = np.frombuffer(f.read(4), dtype=np.uint32)[0]
            spec_win = np.frombuffer(f.read(4), dtype=np.uint32)[0]
            spec_bins = np.frombuffer(f.read(4), dtype=np.uint32)[0]

            rise_time = np.frombuffer(f.read(4), dtype=np.int32)[0]
            peak_amp = np.frombuffer(f.read(4), dtype=np.float32)[0]
            
            amp_envelope = read_float_vector(f)
            centroid = read_float_vector(f)
            flatness = read_float_vector(f)
            
            num_bands = np.frombuffer(f.read(4), dtype=np.uint32)[0]
            band_envelopes_flat = read_float_vector(f)

            # Reshape 1D band envelope data back into 2D (frames x numBands)
            frames = len(band_envelopes_flat) // num_bands if num_bands > 0 else 0
            if frames > 0:
                band_envelopes = band_envelopes_flat.reshape((num_bands, frames)).T
            else: 
                band_envelopes = np.array([])

            transientData = {
                "start": t_start,
                "end": t_end,
                "envHop": env_hop,
                "specHop": spec_hop,
                "specWin": spec_win,
                "specBins": spec_bins,
                "riseTime": rise_time,
                "peakAmp": peak_amp,
                "ampEnvelope": amp_envelope,
                "centroid": centroid,
                "flatness": flatness,
                "numBands": num_bands,
                "bandEnvelopes": band_envelopes
            }
        else:
            # OLD FORMAT: TransientResults (Scalogram)
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

            transientData = {
                "start": t_start,
                "end": t_end,
                "partials": transient_partials
            }

        # Final Block: RMS
        hRms = np.frombuffer(f.read(4), dtype=np.float32)[0]
        tRms = np.frombuffer(f.read(4), dtype=np.float32)[0]

        rmsData = {
            "hRms": hRms,
            "tRms": tRms
        }

    print(f"Loaded {len(indices)} indices, {len(ampData)} fast harmonics.")
    print(f"Transient found: {t_start} to {t_end}.")
    return indices, ampData, freqData, f0, ratios, transientData, rmsData

# =============================
# Plotting
# =============================
def plot_data(indices, fastData, slowData, transientData, 
              plots=["rms", "ratio", "scalogram"], 
              stacked=True,
              amp_threshold=0.0005,
              sr=96000,
              source_separation=False):
    
    if not plots: return

    plt.rcParams.update({
        'font.size': 24, 'axes.titlesize': 24, 'axes.labelsize': 24,
        'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 24,
        'lines.linewidth': 2, 'figure.titlesize': 30, 'font.family': 'sans-serif'
    })

    def to_sec(idx_array):
        return np.array(idx_array) / sr

    axes_list = []
    if stacked:
        num_plots = len(plots)
        fig, axes = plt.subplots(num_plots, 1, figsize=(12, 4 * num_plots), sharex=False)
        if num_plots == 1: axes_list = [axes]
        else: axes_list = list(axes)

    for i, plot_type in enumerate(plots):
        if stacked: ax = axes_list[i]
        else:
            fig = plt.figure(figsize=(10, 5))
            ax = plt.gca()
            fig.canvas.manager.set_window_title(f"Detailed View: {plot_type.upper()}")

        if plot_type == "rms":
            time_axis = to_sec(indices)
            for k, rms in enumerate(fastData):
                plot_len = min(len(time_axis), len(rms))
                ax.plot(time_axis[:plot_len], rms[:plot_len], label=f"H{k+1}")
            ax.set_title("Spectral Envelopes", fontweight='bold')
            ax.set_ylabel("Amplitude")
            ax.set_xlabel("Time (s)")
            ax.grid(True, alpha=0.3)

        elif plot_type in ["freq", "ratio"]:
            data_key = "value" if plot_type == "freq" else "ratio"
            time_axis = to_sec(indices)

            for k, freqPoints in enumerate(slowData):
                if len(freqPoints) == 0: continue
                y_values = np.zeros_like(indices, dtype=float)
                for j, point in enumerate(freqPoints):
                    idx = point["index"]
                    val = point[data_key]
                    is_last_point = (j == len(freqPoints) - 1)
                    next_idx = indices[-1] + 1 if is_last_point else freqPoints[j+1]["index"]
                    mask = (indices >= idx) & (indices < next_idx)
                    y_values[mask] = val
                
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

        elif plot_type == "scalogram":
            t_start = transientData["start"]
            t_end = transientData["end"]
            t_sec_start = t_start / sr
            t_sec_end = t_end / sr
            
            if source_separation:
                # Plot the main amplitude envelope for the new structure
                if "ampEnvelope" in transientData and len(transientData["ampEnvelope"]) > 0:
                    env = transientData["ampEnvelope"]
                    time_vec = np.linspace(t_sec_start, t_sec_end, num=len(env))
                    ax.plot(time_vec, env, color="red", linewidth=2.0)
                    ax.set_title("STransient Amplitude Envelope", fontweight='bold')
            else:
                # Plot partials for the old structure
                for p in transientData["partials"]:
                    if len(p["data"]) == 0: continue
                    time_vec = np.linspace(t_sec_start, t_sec_end, num=len(p["data"]))
                    ax.plot(time_vec, p["data"], alpha=0.7, linewidth=1.5)
                ax.set_title("Transient Partials (Envelope)", fontweight='bold')

            ax.set_ylabel("Amplitude")
            ax.set_xlabel("Time (s)")
            ax.grid(True, alpha=0.3)

        if not stacked: plt.show()
    if stacked: plt.show()

# =============================
# Synthesis Logic
# =============================
def synthesize_transient_old(transientData, sr, total_samples):
    """
    Synthesizes the OLD TransientResults structure (Scalogram sines).
    """
    output = np.zeros(total_samples, dtype=np.float32)
    start = int(transientData["start"])
    end = int(transientData["end"])
    length = end - start
    
    if length <= 0: return output

    t_buffer = np.zeros(length, dtype=np.float32)
    time_indices = np.arange(length)

    for p in transientData["partials"]:
        freq = p["freq"]
        hop = p["hop"]
        env_data = p["data"]
        if len(env_data) == 0: continue

        env_x = np.arange(len(env_data)) * hop
        target_x = np.arange(length)
        interpolator = interp1d(env_x, env_data, kind='linear', bounds_error=False, fill_value=0.0)
        envelope = interpolator(target_x)
        
        phase = 2 * np.pi * freq * time_indices / sr
        sine_wave = np.sin(phase)
        t_buffer += sine_wave * envelope

    write_end = min(end, total_samples)
    write_len = write_end - start
    if write_len > 0:
        output[start:start+write_len] += t_buffer[:write_len]
        
    return output

def synthesize_transient_new(transientData, sr, total_samples):
    output = np.zeros(total_samples, dtype=np.float32)
    start = int(transientData["start"])
    
    env_data = transientData.get("ampEnvelope", [])
    env_hop = transientData.get("envHop", 256)
    
    band_envs = transientData.get("bandEnvelopes", [])
    spec_hop = transientData.get("specHop", 512)
    spec_bins = transientData.get("specBins", 513)
    num_bands = transientData.get("numBands", 0)

    if len(env_data) == 0: 
        return output

    # 1. Determine exact length based on the high-res envelope
    length = len(env_data) * env_hop
    if length <= 0: 
        return output

    # Base Excitation
    base_noise = np.random.normal(0, 1, length).astype(np.float32)
    
    # 2. Spectral Shaping (Filterbank Vocoder)
    if num_bands > 0 and len(band_envs) > 0:
        shaped_noise = np.zeros(length, dtype=np.float32)
        nyquist = sr / 2.0
        bins_per_band = spec_bins / num_bands
        hz_per_bin = nyquist / (spec_bins - 1)
        
        # Time axis for STFT frames
        spec_x = np.arange(len(band_envs)) * spec_hop
        target_x = np.arange(length)
        
        for b in range(num_bands):
            # A. Calculate exact Hz boundaries for this band
            start_bin = b * bins_per_band
            end_bin = spec_bins if (b == num_bands - 1) else (b + 1) * bins_per_band
            
            low_hz = max(20.0, start_bin * hz_per_bin) # Avoid 0Hz DC
            high_hz = min(nyquist - 1.0, end_bin * hz_per_bin)
            
            # B. Bandpass the base noise
            if low_hz < high_hz:
                low_w = low_hz / nyquist
                high_w = high_hz / nyquist
                
                try:
                    # 2nd order Butterworth (gentle crossover)
                    b_coef, a_coef = butter(2, [low_w, high_w], btype='bandpass')
                    # filtfilt ensures zero phase distortion
                    filtered_noise = filtfilt(b_coef, a_coef, base_noise) 
                except ValueError:
                    filtered_noise = np.zeros_like(base_noise)
            else:
                filtered_noise = np.zeros_like(base_noise)
                
            # C. Interpolate the Band Envelope
            band_energy = band_envs[:, b]
            
            # C++ computed mag * mag (Power). We need Amplitude.
            band_amp = np.sqrt(np.maximum(band_energy, 0)) 
            
            interpolator = interp1d(spec_x, band_amp, kind='linear', bounds_error=False, fill_value=0.0)
            smooth_band_env = interpolator(target_x)
            
            # D. Modulate and Accumulate
            shaped_noise += filtered_noise * smooth_band_env
            
    else:
        # Fallback if no band data exists
        shaped_noise = base_noise

    # 3. High-Resolution Temporal Snapping
    # The STFT spectral shaping tracks the "Centroid" and resonance naturally, 
    # but blurs the attack. We apply the temporal envelope as a master VCA.
    env_x = np.arange(len(env_data)) * env_hop
    interpolator = interp1d(env_x, env_data, kind='linear', bounds_error=False, fill_value=0.0)
    master_envelope = interpolator(target_x)
    
    final_transient = shaped_noise * master_envelope

    # 4. Mix into output buffer (allowing it to overlap the harmonics!)
    write_len = min(length, total_samples - start)
    if write_len > 0:
        output[start:start+write_len] += final_transient[:write_len]
        
    return output


def synthesize(f0, ratios, indices, fastData, slowData, transientData, 
               unifiedRMS=None,
               use_transient=True,
               use_harmonics=True,
               source_separation=False,
               sr=48000, filepath="gen/output.wav"):
    
    if not indices.any() and transientData["end"] == 0:
        print("Error: No data found. Cannot synthesize.")
        return

    saveDir = os.path.dirname(filepath)
    if saveDir: os.makedirs(saveDir, exist_ok=True)

    max_harmonic_idx = indices[-1] if len(indices) > 0 else 0
    max_transient_idx = transientData["end"]
    total_samples = max(max_harmonic_idx, max_transient_idx) + 4096 
    
    final_output = np.zeros(total_samples, dtype=np.float32)

    # --- Transient Synthesis ---
    if use_transient:
        print(f"Synthesizing Transient (SourceSep={source_separation})...")
        if source_separation:
            raw_transient = synthesize_transient_new(transientData, sr, total_samples)
        else:
            raw_transient = synthesize_transient_old(transientData, sr, total_samples)

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

    # --- Harmonic Synthesis ---
    if use_harmonics:
        print("Synthesizing Harmonics...")
        harmonic_output = np.zeros(total_samples, dtype=np.float32)
        num_harmonics = min(len(fastData), len(slowData))

        if num_harmonics > 0 and len(indices) > 0:
            for i in range(num_harmonics):
                rms = fastData[i]
                freq_points = slowData[i]
                if len(rms) == 0 or len(freq_points) == 0: continue

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

                freq_segment = np.zeros(range_len, dtype=np.float32)
                curr_freq_val = freq_points[0]["ratio"] * f0
                freq_segment[:] = curr_freq_val 
                
                for k in range(len(freq_points) - 1):
                    p_curr = freq_points[k]
                    p_next = freq_points[k+1]
                    
                    idx_a = p_curr["index"]
                    idx_b = p_next["index"]
                    val_a = p_curr["ratio"] * f0
                    val_b = p_next["ratio"] * f0
                    
                    rel_a = max(0, int(idx_a - start_idx))
                    rel_b = min(range_len, int(idx_b - start_idx))
                    
                    if rel_b > rel_a:
                        ramp = np.linspace(val_a, val_b, rel_b - rel_a, endpoint=False)
                        freq_segment[rel_a : rel_b] = ramp
                
                last_p = freq_points[-1]
                rel_last = int(last_p["index"] - start_idx)
                if rel_last < range_len:
                    freq_segment[rel_last:] = last_p["ratio"] * f0

                phase_segment = np.cumsum(2 * np.pi * freq_segment / sr)
                sine_segment = np.sin(phase_segment)
                harmonic_output[start_idx : start_idx+range_len] += amp_full[start_idx : start_idx+range_len] * sine_segment

        if num_harmonics > 0:
            harmonic_output /= np.sqrt(num_harmonics) 

        if unifiedRMS:
            hStart = int(transientData["end"])
            hEnd = min(len(harmonic_output), hStart + sr)
            h_slice = harmonic_output[hStart:hEnd]
            gen_h_rms = np.sqrt(np.mean(h_slice**2)) if len(h_slice) > 0 else 0.0

            if gen_h_rms > 1e-9 and unifiedRMS["hRms"] > 0:
                gain = unifiedRMS["hRms"] / gen_h_rms
                harmonic_output *= gain
                print(f"   -> Harmonic Gain Applied: {gain:.2f}")

        final_output += harmonic_output
    else:
        print("Skipping Harmonics...")

    # Global Limiter
    peak = np.max(np.abs(final_output))
    if peak > 1.0: final_output = (final_output / peak) * 0.99

    current_peak = np.max(np.abs(final_output))
    target_peak = 0.95

    if current_peak > 1e-9:
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
    instr = "AGuit"
    string = "5"
    freq = "110"
    newF0 = 329.0
    dyn = "F"
    inExt = ".sihat"
    outExt = ".wav"
    outputDir = "sihat/newTests/"
    inputDir = "sihat/newTests/"
    outPrefix = "GEN" + "_"

    plotsStrings=["rms", "ratio", "scalogram"]
    stackplots = False
    ampThreshold = 0.001
    sampleRate = 96000
    replaceF0 = False
    
    if replaceF0: outPrefix = "GENto" + str(int(newF0)) + "_"
    
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
    parser.add_argument("--source_sep", action="store_true", default=True, help="Enable parsing of STransientResults structure")
    
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
                    indices, fast, slow, f0, ratios, tData, rmsData = load_sihat_data(file_path, args.source_sep)
                    
                    rawName = file_path.stem + ".wav"
                    output_filename = rawName[len(inPrefix):] if rawName.startswith(inPrefix) else rawName
                    target_f0 = newF0 if replaceF0 else f0
                    
                    synthesize(target_f0, ratios, indices, fast, slow, tData, 
                               use_harmonics=useOvertones, use_transient=useTransient, 
                               source_separation=args.source_sep,
                               sr=args.sr, filepath=str(output_dir / output_filename), unifiedRMS=rmsData)
                except Exception as e:
                    print(f"Failed: {e}")

        elif args.mode in ["plot", "synth"]:
            if not input_path.is_file():
                print(f"Error: input_path must be a file.")
                sys.exit(1)
                
            print(f"Reading: {input_path}")
            indices, fast, slow, f0, ratios, tData, rmsData = load_sihat_data(input_path, args.source_sep)
            
            if args.mode == "plot":
                plot_data(indices, fast, slow, tData, plotsStrings, stackplots, ampThreshold, source_separation=args.source_sep)
            else:
                target_f0 = newF0 if replaceF0 else f0
                synthesize(target_f0, ratios, indices, fast, slow, tData, 
                           use_harmonics=useOvertones, use_transient=useTransient, 
                           source_separation=args.source_sep,
                           sr=args.sr, filepath=str(args.out), unifiedRMS=rmsData)
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)