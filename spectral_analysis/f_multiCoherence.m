function [Cxy,phi,f] = f_multiCoherence(sig,fs,varargin)

p = inputParser;
addParameter(p,'f_res',[]);
addParameter(p,'K',15);

parse(p,varargin{:});

% Adjust segment length based on sampling frequency
alpha = 0.5;
K = p.Results.K;

T = size(sig,1);
segment_length = 2^nextpow2(T / (0.5*K + 0.5));

% Parameters for Welch's method
overlap = alpha*segment_length;        % Overlap between segments
nfft = segment_length;                        % Number of FFT points

% Divide signals into overlapping segments
window = hann(segment_length);
step = segment_length - overlap;
num_segments = floor((size(sig,1) - overlap) / step);

N = size(sig,2);

% Initialize accumulators for spectral densities
Pxx_accum = zeros(nfft, N);
Pxy_accum = zeros(nfft, N, N);
%%
for i = 1:num_segments
    start_idx = (i-1) * step + 1;
    end_idx = start_idx + segment_length - 1;

    % Extract segments and apply window
    segment1 = sig(start_idx:end_idx,:) .* window;
    segment2 = sig(start_idx:end_idx,:) .* window;

    % Compute FFTs
    fft_segment = fft(segment1, nfft);

    % Accumulate power spectral densities and cross-spectral density
    Pxx_accum = Pxx_accum + (1/(fs*segment_length)) * abs(fft_segment).^2;
    Pxy_accum = Pxy_accum + (1/(fs*segment_length)) * (fft_segment .* conj(permute(fft_segment,[1,3,2])));
end

%%
% Average the spectral densities
Pxx = Pxx_accum / num_segments;
Pxy = Pxy_accum / num_segments;

% Compute the coherence
Cxy = abs(Pxy).^2 ./ (Pxx .* permute(Pxx,[1,3,2]));

% Compute the phase
phi = angle(Pxy);

% Frequency vector
f = (0:nfft-1) * (fs/nfft);

% Only keep the positive frequencies
half_nfft = floor(nfft / 2) + 1;
Cxy = permute(Cxy(1:half_nfft,:,:),[2,3,1]);
f = f(1:half_nfft);

% Only keep the positive frequencies for the phase
phi = permute(phi(1:half_nfft,:,:), [2,3,1]);

end