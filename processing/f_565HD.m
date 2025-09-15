function [rfp,rfp_HD] = f_565HD(rfp,HbO,HbR)

rfp = rfp./mean(rfp,3);

if nargin > 1
    tmpHbRed = HbR./mean(HbR,3);
    tmpHbGreen = HbO./mean(HbO,3);
    rfp_HD = rfp./(tmpHbRed.^0.8.*tmpHbGreen.^0.4);
    rfp_HD = rfp_HD-1;
else
    rfp_HD = [];
end

rfp = rfp-1;
