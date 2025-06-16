import librosa
import numpy as np
import soundfile as sf

def split_audio_into_periods(input_file, output_prefix="guit"):
    # Load audio file
    y, sr = librosa.load(input_file, sr=None, mono=True)
    
    # Compute amplitude envelope
    envelope = np.abs(y)
    window_size = int(0.01 * sr)  # 10ms window for smoothing
    envelope = np.convolve(envelope, np.ones(window_size)/window_size, mode='same')
    
    # Detect sustain phase
    max_amp = np.max(envelope)
    threshold_onset = 0.2 * max_amp
    threshold_sustain = 0.2 * max_amp
    sustain_start = None
    sustain_end = None

    # Find sustain start
    for i in range(len(envelope)):
        if envelope[i] >= threshold_onset:
            N = int(0.05 * sr)  # 50ms
            if i + N >= len(envelope):
                break
            if np.all(envelope[i:i+N] >= threshold_sustain):
                sustain_start = i
                break

    # Find sustain end
    for i in reversed(range(len(envelope))):
        if envelope[i] >= threshold_sustain:
            N = int(0.05 * sr)
            start = max(0, i - N)
            if np.all(envelope[start:i] >= threshold_sustain):
                sustain_end = i
                break

    if sustain_start is None or sustain_end is None:
        raise ValueError("Could not detect sustain portion in the audio.")
    
    sustain_signal = y[sustain_start:sustain_end]
    
    # Estimate fundamental frequency using YIN on a stable segment
    segment_start = len(sustain_signal) // 2
    segment_length = 2048
    if segment_start + segment_length > len(sustain_signal):
        segment_length = len(sustain_signal) - segment_start
    segment = sustain_signal[segment_start:segment_start + segment_length]
    
    f0 = librosa.yin(segment, fmin=50, fmax=2000, sr=sr)
    f0 = f0[f0 > 0]  # Filter out invalid values
    if len(f0) == 0:
        raise ValueError("Fundamental frequency detection failed.")
    average_f0 = np.median(f0)
    period_samples = int(sr / average_f0)
    
    # Ensure period_samples is valid
    if period_samples <= 0 or period_samples > len(sustain_signal) // 2:
        raise ValueError("Invalid period detected.")
    
    # Truncate to integer number of periods
    num_periods = len(sustain_signal) // period_samples
    truncated_length = num_periods * period_samples
    sustain_truncated = sustain_signal[:truncated_length]
    
    # Split into periods
    periods = np.array_split(sustain_truncated, num_periods)
    
    # Save each period
    for i, period in enumerate(periods):
        output_file = f"{output_prefix}_{i}.wav"
        sf.write("./media/" + output_file, period, sr)
    
    return len(periods)

# Example usage
if __name__ == "__main__":
    input_audio = "./media/EGuit_Low_E.wav"
    try:
        num_periods = split_audio_into_periods(input_audio)
        print(f"Successfully split into {num_periods} periods.")
    except Exception as e:
        print(f"Error: {e}")