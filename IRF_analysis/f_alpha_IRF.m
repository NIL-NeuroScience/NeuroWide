%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            f_alpha_IRF
% author - Brad Rauscher (created 2024)
% 
% Creates IRF using a double gamma function.
% 
% INPUTS: f_alpha_IRF(t0, tau1, tau2, A, B, range)
%   t0: temporal offset
%   tau1: decay constant for first gamma function
%   tau2: decay constant for second gamma function
%   A: amplitude for first gamma function
%   B: amplitude for second gamma function
%   range: time window [t1, t2]
% 
% OUTPUTS:
%   IRF: resulting double gamma function
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IRF = f_alpha_IRF(t0, tau1, tau2, A, B, range)

    tr = ((range(1) : range(2)) - t0)';
    D = (tr ./ tau1).^3 .* exp(-tr ./ tau1);
    D(tr < 0) = 0;
    C = (tr ./ tau2).^3 .* exp(-tr ./ tau2);
    C(tr < 0) = 0;
    
    IRF = A .* D + B .* C;
end