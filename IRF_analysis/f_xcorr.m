function [c] = f_xcorr(sig1,sig2,maxlag)
% calculates cross-correlation function between columns of sig1 and sig2
% 
% Inputs:
%   sig1 - 1st signal (t x N1)
%   sig2 - 2nd signal (t x N2)
%   maxlag - max negative and positive lag (datapoints)

%%
% sig1 = Signals.Mid(2,:)';
% sig2 = Signals.Mid(5,:)';
% maxlag = 30;
% end of parameters

sig1 = sig1-mean(sig1,1);
sig2 = sig2-mean(sig2,1);

m = numel(sig1(:,1));
maxlagDefault = m-1;
mxl = min(maxlag,maxlagDefault);
m2 = findTransformLength(m);

% tapers = [3,1];
% tapers = dpss(m,tapers(1),tapers(2));
% tapers = tapers*sqrt(10);
% tapers = permute(tapers,[1 3 2]);

% X = sig1.*tapers;
% Y = sig2.*tapers;
X = sig1;
Y = sig2;

X = fft(X,m2,1);
Y = fft(Y,m2,1);

% X = mean(X,3);
% Y = mean(Y,3);

c1 = ifft(X.*conj(Y),[],1,'symmetric');

c = [c1(m2 - mxl + (1:mxl),:); c1(1:mxl+1,:)];

% X = mean(sig1.*tapers,3);
% Y = mean(sig2.*tapers,3);
X = sig1;
Y = sig2;
cxx0 = sum(X.^2,1);
cyy0 = sum(Y.^2,1);
scaleCoeffCross = sqrt(cxx0.*cyy0);

c = c./scaleCoeffCross;

% c = mean(c,3);
%%
%-------------------------------------------------------------------------
function m = findTransformLength(m)
m = 2*m;
while true
    r = m;
    for p = [2 3 5 7]
        while (r > 1) && (mod(r, p) == 0)
            r = r / p;
        end
    end
    if r == 1
        break;
    end
    m = m + 1;
end
