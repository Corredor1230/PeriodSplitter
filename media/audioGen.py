import numpy as np
import pandas as pd
from scipy.io.wavfile import write

def synthesize_simple(amps_df, samples, Fs):
    """
    amps_df: DataFrame with frequencies as headers, amplitudes as rows
    samples: array of sample indices corresponding to each row in amps_df
    Fs: sample rate
    """
    freqs = amps_df.columns.astype(float).to_numpy()   # partial frequencies
    amps = amps_df.to_numpy()                          # shape (frames, partials)
    B, P = amps.shape

    N = samples[-1]  # assume last sample is end
    t = np.arange(N) / Fs
    x = np.zeros(N, dtype=np.float64)

    # For each partial, generate continuous sine and apply amplitude envelope
    for k, f in enumerate(freqs):
        sine = np.cos(2 * np.pi * f * t)

        # build amplitude envelope by interpolating between checkpoints
        env = np.zeros(N)
        for b in range(len(samples)-2 or len(amps) -2):
            i0, i1 = samples[b], samples[b+1]
            if i1 <= i0:
                continue
            a0, a1 = amps[b, k], amps[b+1, k]
            # linear interpolation of amplitude
            env[i0:i1] = np.linspace(a0, a1, i1-i0)

        x += sine * env

    return x

# ==== Example usage ====
Fs = 96000
fileName = "EGuit6_82_F"

amps_df = pd.read_csv("AMP" + fileName + ".csv")
samples = pd.read_csv(fileName + ".csv", header=None).values.squeeze()

x = synthesize_simple(amps_df, samples, Fs)

# Normalize and save
# x /= np.max(np.abs(x)) + 1e-12
write("REC" + fileName + ".wav", Fs, (x*32767).astype(np.int16))

print("Saved")
