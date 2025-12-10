function IRF = f_alpha_IRF(t0,tau1,tau2,A,B,range)
    hrf_l = range;
    tr = (hrf_l-t0)';
    D = ((tr)./tau1).^3.*exp(-(tr)./tau1);
    D(tr<0) = 0;
    C = ((tr)./tau2).^3.*exp(-(tr)./tau2);
    C(tr<0) = 0;
    
    IRF = A.*D + B.*C;
end

    % hrf_l = range;
    % tr = (((hrf_l(1):hrf_l(2)))-t0)';
    % D = ((tr)./tau1).^3.*exp(-(tr)./tau1);
    % D(tr<0) = 0;
    % C = ((tr)./tau2).^3.*exp(-(tr)./tau2);
    % C(tr<0) = 0;
    % 
    % IRF = A.*D + B.*C;