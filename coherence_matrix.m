% Function to calculate coherence between two time signals
function [Cxy, f] = calculate_coherence(signal1, signal2, fs)
    % Inputs:
    % signal1 - First time signal
    % signal2 - Second time signal
    % fs - Sampling frequency (Hz)
    %
    % Outputs:
    % Cxy - Coherence values
    % f - Frequency vector

    % Check if the signals are of the same length
    if length(signal1) ~= length(signal2)
        error('Signals must have the same length.');
    end

    % Define parameters for the coherence calculation
    window = hamming(256); % Hamming window of length 256
    noverlap = 128;        % 50% overlap
    nfft = 512;            % Number of FFT points

    % Calculate the coherence
    [Cxy, f] = mscohere(signal1, signal2, window, noverlap, nfft, fs);

    % Plot the coherence
    figure;
    plot(f, Cxy);
    xlabel('Frequency (Hz)');
    ylabel('Coherence');
    title('Coherence between Signal 1 and Signal 2');
    grid on;
end

% Function to calculate coherence between two time signals from scratch
function [Cxy, f] = calculate_coherence_from_scratch(signal1, signal2, fs)
    % Inputs:
    % signal1 - First time signal
    % signal2 - Second time signal
    % fs - Sampling frequency (Hz)
    %
    % Outputs:
    % Cxy - Coherence values
    % f - Frequency vector

    % Check if the signals are of the same length
    if length(signal1) ~= length(signal2)
        error('Signals must have the same length.');
    end

    % Length of the signals
    N = length(signal1);

    % Compute the Fourier transforms of the signals
    fft_signal1 = fft(signal1);
    fft_signal2 = fft(signal2);

    % Compute the power spectral densities
    Pxx = (1/(fs*N)) * abs(fft_signal1).^2; % Power spectral density of signal1
    Pyy = (1/(fs*N)) * abs(fft_signal2).^2; % Power spectral density of signal2

    % Compute the cross-spectral density
    Pxy = (1/(fs*N)) * (fft_signal1 .* conj(fft_signal2));

    % Compute the coherence
    Cxy = abs(Pxy).^2 ./ (Pxx .* Pyy);

    % Frequency vector
    f = (0:N-1) * (fs/N);

    % Plot the coherence
    figure;
    plot(f, Cxy);
    xlabel('Frequency (Hz)');
    ylabel('Coherence');
    title('Coherence between Signal 1 and Signal 2 (from scratch)');
    grid on;
end