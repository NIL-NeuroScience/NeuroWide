function [Hb,HbO] = f_calcHb(ch525,ch625,varargin)

p = inputParser;
addParameter(p,'dim',3);
addParameter(p,'ch525_baseline',[]);
addParameter(p,'ch625_baseline',[]);

parse(p,varargin{:});

if ~isempty(p.Results.ch525_baseline)
    ch525_baseline = mean(p.Results.ch525_baseline,p.Results.dim);
else
    ch525_baseline = mean(ch525,p.Results.dim);
end

if ~isempty(p.Results.ch625_baseline)
    ch625_baseline = mean(p.Results.ch625_baseline,p.Results.dim);
else
    ch625_baseline = mean(ch625,p.Results.dim);
end

settings.hd=f_defineHemodynamicParameters(525,625);

tmpA0Hb=settings.hd.cLambda2Hb*log(ch625_baseline)-settings.hd.cLambda1Hb*log(ch525_baseline);
tmpA0HbO=settings.hd.cLambda2HbO*log(ch625_baseline)-settings.hd.cLambda1HbO*log(ch525_baseline);
HbO=(tmpA0HbO+settings.hd.cLambda1HbO*log(double(ch525))-settings.hd.cLambda2HbO*log(double(ch625)));
Hb=(tmpA0Hb+settings.hd.cLambda1Hb*log(double(ch525))-settings.hd.cLambda2Hb*log(double(ch625)));
