import numpy as np
import pandas as pd
from scipy.io.wavfile import write

def synthesize_simple(amps_df, Fs):
    """
    amps_df: DataFrame with frequencies as headers, amplitudes as rows,
             and sample indices as the index
    Fs: sample rate
    """
    freqs = amps_df.columns.astype(float).to_numpy()   # partial frequencies
    amps = amps_df.to_numpy()                          # shape (frames, partials)
    samples = amps_df.index.to_numpy()                 # sample indices
    B, P = amps.shape

    N = samples[-1]  # assume last sample index is end
    t = np.arange(N) / Fs
    x = np.zeros(N, dtype=np.float64)

    # For each partial, generate continuous sine and apply amplitude envelope
    for k, f in enumerate(freqs):
        sine = np.cos(2 * np.pi * f * t)

        env = np.zeros(N)
        for b in range(len(samples)-1):
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
fileName = "AGuit4_146_F"

# Load amps_df with samples as row headers
amps_df = pd.read_csv("AMP" + fileName + ".csv", index_col=0)

x = synthesize_simple(amps_df, Fs)

# Save
from scipy.io.wavfile import write
write("REC" + fileName + ".wav", Fs, (x/np.max(np.abs(x)) * 32767).astype(np.int16))

print("Saved")
