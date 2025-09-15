function [Hb,HbO] = f_calcHb(ch525,ch625)

settings.hd=f_defineHemodynamicParameters_WF(525,625);

tmpA0Hb=settings.hd.cLambda2Hb*log(mean(ch625,3))-settings.hd.cLambda1Hb*log(mean(ch525,3));
tmpA0HbO=settings.hd.cLambda2HbO*log(mean(ch625,3))-settings.hd.cLambda1HbO*log(mean(ch525,3));
HbO=(tmpA0HbO+settings.hd.cLambda1HbO*log(double(ch525))-settings.hd.cLambda2HbO*log(double(ch625)));
Hb=(tmpA0Hb+settings.hd.cLambda1Hb*log(double(ch525))-settings.hd.cLambda2Hb*log(double(ch625)));

