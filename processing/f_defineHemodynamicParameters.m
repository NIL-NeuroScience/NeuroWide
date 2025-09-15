function hd = f_defineHemodynamicParameters(Lambda1,Lambda2)
hd.Lambda1=Lambda1;
hd.Lambda2=Lambda2;
tmp1=f_GetExtinctions_WF([hd.Lambda1 hd.Lambda2]);
tmp2=f_pathlengths_WF([hd.Lambda1 hd.Lambda2]);
hd.exsLambda1Hb=tmp1(1,2);
hd.exsLambda1HbO=tmp1(1,1);
hd.exsLambda2Hb=tmp1(2,2);
hd.exsLambda2HbO=tmp1(2,1);
hd.pathLambda1=tmp2(1);
hd.pathLambda2=tmp2(2);
hd.cLambda1HbO=1/hd.pathLambda1*(hd.exsLambda2Hb/(hd.exsLambda2HbO*hd.exsLambda1Hb-hd.exsLambda2Hb*hd.exsLambda1HbO));
hd.cLambda2HbO=1/hd.pathLambda2*(hd.exsLambda1Hb/(hd.exsLambda2HbO*hd.exsLambda1Hb-hd.exsLambda2Hb*hd.exsLambda1HbO));
hd.cLambda1Hb=1/hd.pathLambda1*(hd.exsLambda2HbO/(hd.exsLambda2Hb*hd.exsLambda1HbO-hd.exsLambda2HbO*hd.exsLambda1Hb));
hd.cLambda2Hb=1/hd.pathLambda2*(hd.exsLambda1HbO/(hd.exsLambda2Hb*hd.exsLambda1HbO-hd.exsLambda2HbO*hd.exsLambda1Hb));
end