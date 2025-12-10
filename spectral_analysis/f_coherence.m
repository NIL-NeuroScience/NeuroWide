% Function to calculate coherence between two time signals using Welch's method
function [Cxy, f] = f_coherence(signal1, signal2, fs, varargin)
    % Inputs:
    % signal1 - First time signal
    % signal2 - Second time signal
    % fs - Sampling frequency (Hz)
    %
    % Outputs:
    % Cxy - Coherence values
    % f - Frequency vector

    p = inputParser;
    addParameter(p,'method','welch'); % or 'chronux'
    addParameter(p,'tapers',[5,9]);
    
    parse(p,varargin{:});

    % Check if the signals are of the same length
    if length(signal1) ~= length(signal2)
        error('Signals must have the same length.');
    end
    
    if string(p.Results.method) == "welch"
    
        % Adjust segment length based on sampling frequency
        base_segment_length = 256; % Base segment length for a standard sampling frequency
        segment_length = max(base_segment_length, round(fs / 10)); % Adjust segment length to be proportional to fs
    
        % Ensure segment length is a power of 2 for FFT efficiency
        segment_length = 2^nextpow2(segment_length);
    
        % Parameters for Welch's method
        overlap = 128;        % Overlap between segments
        nfft = 512;           % Number of FFT points
    
        % Divide signals into overlapping segments
        window = hamming(segment_length);
        step = segment_length - overlap;
        num_segments = floor((length(signal1) - overlap) / step);
    
        % Initialize accumulators for spectral densities
        Pxx_accum = zeros(nfft, 1);
        Pyy_accum = zeros(nfft, 1);
        Pxy_accum = zeros(nfft, 1);
    
        % Loop over segments
        for i = 1:num_segments
            start_idx = (i-1) * step + 1;
            end_idx = start_idx + segment_length - 1;
    
            % Extract segments and apply window
            segment1 = signal1(start_idx:end_idx) .* window;
            segment2 = signal2(start_idx:end_idx) .* window;
    
            % Compute FFTs
            fft_segment1 = fft(segment1, nfft);
            fft_segment2 = fft(segment2, nfft);
    
            % Accumulate power spectral densities and cross-spectral density
            Pxx_accum = Pxx_accum + (1/(fs*segment_length)) * abs(fft_segment1).^2;
            Pyy_accum = Pyy_accum + (1/(fs*segment_length)) * abs(fft_segment2).^2;
            Pxy_accum = Pxy_accum + (1/(fs*segment_length)) * (fft_segment1 .* conj(fft_segment2));
        end
    
        % Average the spectral densities
        Pxx = Pxx_accum / num_segments;
        Pyy = Pyy_accum / num_segments;
        Pxy = Pxy_accum / num_segments;
    
        % Compute the coherence
        Cxy = abs(Pxy).^2 ./ (Pxx .* Pyy);
    
        % Frequency vector
        f = (0:nfft-1) * (fs/nfft);
    
        % Only keep the positive frequencies
        half_nfft = floor(nfft / 2) + 1;
        Cxy = Cxy(1:half_nfft);
        f = f(1:half_nfft);

    else
        params = struct;
        params.Fs = fs;
        params.tapers = p.Results.tapers;
        params.trialave = 0;
        
        [Cxy,phi,~,~,~,f] = coherencyc(signal1,signal2,params);
        
        f = f';
    end
end