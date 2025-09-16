function PathLengths = f_pathlengths(lambda,factor)
% in cm
if nargin == 1
    factor = 1;
end
Hb0 = 60*10^-6;
Hb = 40*10^-6;
% lambda = [473, 530, 590, 625];
g = 0.9;
c = 3*10^10;
e = f_GetExtinctions(lambda);
mua = e(:,1).*Hb0+e(:,2).*Hb;
mus = 150*(lambda/560).^(-2);
mus = mus*factor;
z0 = 1./((1-g)*mus');
gamma = sqrt(c./(3*(mua+(1-g)*mus')));
PathLengths = (c*z0./(2*gamma.*(mua*c).^0.5)).*(1+(3/c)*mua.*(gamma).^2);