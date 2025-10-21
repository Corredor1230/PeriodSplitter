import argparse
import numpy as np
import matplotlib.pyplot as plt
import soundfile as sf
import os

# =============================
# Binary Loaders
# =============================
def load_indices(path):
    with open(path, "rb") as f:
        n = np.frombuffer(f.read(4), dtype=np.uint32)[0]
        data = np.frombuffer(f.read(n * 4), dtype=np.uint32)
    return data

def load_fast(path):
    with open(path, "rb") as f:
        num_harmonics = np.frombuffer(f.read(4), dtype=np.uint32)[0]
        harmonics = []
        for _ in range(num_harmonics):
            length = np.frombuffer(f.read(4), dtype=np.uint32)[0]
            data = np.frombuffer(f.read(length * 4), dtype=np.float32)
            harmonics.append(data)
    return harmonics

def load_slow(path):
    with open(path, "rb") as f:
        num_harmonics = np.frombuffer(f.read(4), dtype=np.uint32)[0]
        harmonics = []
        for _ in range(num_harmonics):
            length = np.frombuffer(f.read(4), dtype=np.uint32)[0]
            entries = np.frombuffer(f.read(length * 12), dtype=[("index", np.uint32), ("value", np.float32), ("ratio", np.float32)])
            harmonics.append(entries)
    return harmonics

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
        axes[0].plot(indices[:len(rms)], rms, label=f"H{i+1}")
    axes[0].set_title("Fast RMS (per harmonic)")
    axes[0].set_ylabel("Amplitude (RMS)")
    axes[0].legend(loc="upper right", ncol=4, fontsize=8)
    axes[0].grid(True, alpha=0.3)

    # ---------------------------
    # 2️⃣ Frequency per harmonic
    # ---------------------------
    for i, freqPoints in enumerate(slowData):
        freq = np.zeros_like(indices, dtype=float)
        for j, (idx, val, ratio) in enumerate(zip(freqPoints["index"], freqPoints["value"], freqPoints["ratio"])):
            next_idx = indices[-1] if j == len(freqPoints)-1 else freqPoints["index"][j+1]
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
        ratio = np.zeros_like(indices, dtype=float)
        for j, (idx, val, r) in enumerate(zip(freqPoints["index"], freqPoints["value"], freqPoints["ratio"])):
            next_idx = indices[-1] if j == len(freqPoints)-1 else freqPoints["index"][j+1]
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
# Synth mode
# =============================
def synthesize(indices, fastData, slowData, sr=48000, filename="default", uniqueDirectory="gen/"):
    # Linear interpolation helpers
    def interp(vec, size):
        return np.interp(np.arange(size), np.linspace(0, size, len(vec)), vec)

    savePath = uniqueDirectory + filename + ".wav"
    total_samples = indices[-1]
    t = np.arange(total_samples) / sr
    output = np.zeros_like(t)

    for i, (rms, freq_points) in enumerate(zip(fastData, slowData)):
        if len(rms) == 0 or len(freq_points) == 0: continue
        amp = np.interp(np.arange(len(indices)), np.linspace(0, len(indices)-1, len(rms)), rms)
        freq = np.zeros_like(indices, dtype=float)
        for j, (idx, val, r) in enumerate(freq_points):
            next_idx = indices[-1] if j == len(freq_points)-1 else freq_points[j+1]["index"]
            mask = (indices >= idx) & (indices < next_idx)
            start_val = val
            end_val = val if j == len(freq_points)-1 else freq_points[j+1]["value"]
            interp_freq = np.linspace(start_val, end_val, mask.sum(), endpoint=False)
            freq[mask] = interp_freq

        # Interpolate freq and amp to match total sample length
        freq_full = np.interp(np.arange(total_samples), indices, freq)
        amp_full = np.interp(np.arange(total_samples), indices, amp)

        phase = np.cumsum(2 * np.pi * freq_full / sr)
        output += amp_full * np.sin(phase)

    output /= len(fastData)
    # sd.play(output, sr)
    # sd.wait()
    os.makedirs(os.path.dirname(savePath), exist_ok=True)

    # Normalization
    target_peak = 0.9
    peak = np.max(np.abs(output))
    if peak > 0:
        output = (output / peak) * target_peak

    sf.write(savePath, output, sr)

# =============================
# Main
# =============================
if __name__ == "__main__":

    instr = "EGuit6"
    freq = "82"
    dyn = "F"
    flag = "plot"

    filename = instr + "_" + freq + "_" + dyn
    outname = "GEN_" + filename
    indName = "INDEX" + filename + ".bin"
    fasName = "AMP" + filename + ".bin"
    sloName = "FREQ" + filename + ".bin"

    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["plot", "synth"], default=flag)
    parser.add_argument("--index", default=indName)
    parser.add_argument("--fast", default=fasName)
    parser.add_argument("--slow", default=sloName)
    args = parser.parse_args()

    indices = load_indices(args.index)
    fastData = load_fast(args.fast)
    slowData = load_slow(args.slow)

    if args.mode == "plot":
        plot_data(indices, fastData, slowData)
    else:
        synthesize(indices, fastData, slowData, 96000, filename=outname)
