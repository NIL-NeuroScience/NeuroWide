function [spg,t,fSpectogram] = f_morlet(y,fs,fwin,fsteps)
%%

t = 1/fs:1/fs:size(y,1)/fs;

fSpectogram=logspace(log10(fwin(1)),log10(fwin(2)),fsteps);
morletFWHM=5;

%%
Ts = 1/fs;

sigma_tc = morletFWHM / sqrt(8*log(2));
sigma_t = sigma_tc ./ fSpectogram;
nscales = length(fSpectogram);
precision = 3;

nx = size(y,1);
ntimes = size(y,2);
P = zeros(nx,nscales,ntimes);

for s = 1:nscales
    xval = (-precision*sigma_t(s) : 1/fs : precision*sigma_t(s))';
    W = sqrt(fSpectogram(s)) * morlet_wavelet(fSpectogram(s)*xval, sigma_tc);
    P(:,s,:) = conv2(y, W, 'same') * Ts; 
end

spg = abs(P).^2;

%%

function W = morlet_wavelet(t,sigma_tc)
    W = (sigma_tc*sqrt(pi))^(-0.5) * exp( -(t.^2)/(2*sigma_tc^2) ) .* exp(1i*2*pi*t);
end

end