import soundfile as sf
import numpy as np
from numpy.fft import fft, ifft
import matplotlib.pyplot as plt  # These are the imports which we needed
import librosa
from scipy import signal


def pitch_shift(signal, fs, f_ratio):
    '''
    This function helps in shifting the pitch in the time domain by a factor "f_ratio" 
    f_ratio is ratio between new frequncy needed and original frequency
    '''
    peaks = find_peaks(signal, fs)
    new_signal = psola(signal, peaks, f_ratio)
    return new_signal


def find_peaks(signal, fs, max_hz=950, min_hz=75, analysis_win_ms=40, max_change=1.05, min_change=0.95):
    '''
    The function first computes the pitch periodicity by dividing the signal into small overlapping windows and 
    computing the period of each window.Then it finds the peaks in the signal that correspond to the fundamental 
    frequency by searching for the maximum values near the expected location. The function returns an array of the peak 
    positions in the input signal.
    '''
    N = len(signal)
    min_period = fs // max_hz
    max_period = fs // min_hz

    
    sequence = int(analysis_win_ms / 1000 * fs)
    periods = periods_per_seq(signal, sequence, min_period, max_period)

    peaks = [np.argmax(signal[:int(periods[0]*1.1)])]
    while True:
        prev = peaks[-1]
        idx = prev // sequence 
        if prev + int(periods[idx] * max_change) >= N:
            break
        peaks.append(prev + int(periods[idx] * min_change) + np.argmax(signal[prev + int(periods[idx] * min_change): prev + int(periods[idx] * max_change)]))
    return np.array(peaks)


def periods_per_seq(signal, sequence, min_period, max_period):
    '''
    The function computes the pitch periodicity of the signal by dividing it into small windows 
    and computing the period of each window using autocorrelation.The autocorrelation peak period is found by computing
    the FFT of each window, setting the DC component to zero, taking the inverse FFT, and finding the maximum value in 
    the specified period range.The function returns an array of the autocorrelation peak periods for each window sequence.
    '''
    offset = 0  # current sample offset
    periods = []  # period length of each analysis sequence

    while offset < N:
        fourier = fft(signal[offset: offset + sequence])
        fourier[0] = 0  # remove DC component
        auto_correlation = ifft(fourier * np.conj(fourier)).real
        auto_correlation_peak = min_period + \
            np.argmax(auto_correlation[min_period: max_period])
        periods.append(auto_correlation_peak)
        offset += sequence
    return periods


def psola(signal, peaks, f_ratio):
    '''
     The function first performs linear interpolation to change the number of peak indices to match the desired 
     time-stretch factor. Then by iterating over each new peak index, finding the corresponding old peak index, and 
     constructing an overlapping window centered at the new peak that is applied to a portion of the input signal 
     centered at the old peak. The resulting windowed signals are then added together to produce the output signal with 
     the desired time-stretch factor. This functions return the final sginal.
    '''
    N = len(signal)

    new_signal = np.zeros(N)
    print("The value of f_ratio is", f_ratio)

    new_peaks_ref = np.linspace(0, len(peaks) - 1, int(len(peaks) * f_ratio))
    new_peaks = np.zeros(len(new_peaks_ref)).astype(int)  

    for i in range(len(new_peaks)):
        weight = new_peaks_ref[i] % 1
        left = np.floor(new_peaks_ref[i]).astype(int)
        right = np.ceil(new_peaks_ref[i]).astype(int)
        new_peaks[i] = int(peaks[left] * (1 - weight) + peaks[right] * weight)


    for j in range(len(new_peaks)):
        # find the corresponding old peak index
        i = np.argmin(np.abs(peaks - new_peaks[j]))
        # get the distances to adjacent peaks
        P1 = [new_peaks[j] if j == 0 else new_peaks[j] - new_peaks[j-1],N - 1 - new_peaks[j] if j == len(new_peaks) - 1 else new_peaks[j+1] - new_peaks[j]]
        # edge case truncation
        if peaks[i] - P1[0] < 0:
            P1[0] = peaks[i]
        if peaks[i] + P1[1] > N - 1:
            P1[1] = N - 1 - peaks[i]
        # linear OLA window
        window = list(np.linspace(0, 1, P1[0] + 1)[1:]) + list(np.linspace(1, 0, P1[1] + 1)[1:])
        # center window from original signal at the new peak
        new_signal[new_peaks[j] - P1[0]: new_peaks[j] + P1[1]] += window * signal[peaks[i] - P1[0]: peaks[i] + P1[1]]
    return new_signal






sample_rate = 44100                            # These are parameters for lowpass filter
cutoff_freq = 2000  # Hz
nyquist_freq = 0.5 * sample_rate
cutoff_normalized = cutoff_freq / nyquist_freq
order = 5


b, a = signal.butter(order, cutoff_normalized, btype='low')   # Create the Butterworth filter

orig_signal, fs = librosa.load("Voice2.wav", sr=44100)            #loading the voice

orig_signal = signal.lfilter(b, a, orig_signal)                    #using the LPF filter
N = len(orig_signal)

f_ratio = 0.6                                                       # Pitch shift amount as a ratio

new_signal = pitch_shift(orig_signal, fs, f_ratio)                  # Shift pitch

new_signal = 10*new_signal                                          #Increase the amplitude
 
sf.write('female_scale_transposed.wav', new_signal, fs)              #Write the output
