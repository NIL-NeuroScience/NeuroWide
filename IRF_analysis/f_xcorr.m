%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               f_xcorr
% author - Brad Rauscher (created 2024)
% 
% Calculates cross-correlation between columns of matrices 'sig1' and 
% 'sig2'. Removes mean of each column in 'sig1' and 'sig2'.
% 
% INPUTS: f_xcorr(sig1, sig2, maxlag)
%   sig1: first signal (T x N)
%   sig2: second signal (T x N)
%   maxlag: max negative and positive lag
% 
% OUTPUTS:
%   r: correlation value(s) for each lag
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = f_xcorr(sig1, sig2, maxlag)

% remove mean
sig1 = sig1 - mean(sig1, 1);
sig2 = sig2 - mean(sig2, 1);

% calculate maxlag and find transform length
N = size(sig1, 1);
maxlagDefault = N - 1;
mxl = min(maxlag, maxlagDefault);
m2 = findTransformLength(N);

% calculate cross-correlation
X = fft(sig1, m2, 1);
Y = fft(sig2, m2, 1);

c1 = ifft(X .* conj(Y), [], 1, 'symmetric');

% correctly index cross-correlation
r = [c1(m2 - mxl + (1 : mxl), :); c1(1 : mxl + 1, :)];

% rescale to pearson's coefficient
cxx0 = sum(sig1.^2, 1);
cyy0 = sum(sig2.^2, 1);
scaleCoeffCross = sqrt(cxx0 .* cyy0);

r = r ./ scaleCoeffCross;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRA FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m = findTransformLength(m)

m = 2 * m;

while true
    r = m;
    for p = [2, 3, 5, 7]
        while (r > 1) && (mod(r, p) == 0)
            r = r / p;
        end
    end
    if r == 1
        break;
    end
    m = m + 1;
end
