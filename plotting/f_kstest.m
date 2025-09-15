function [h,p] = f_kstest(X,alpha)

N = numel(X);

h = zeros(N);
p = zeros(N);

for hI = 1:N
    for wI = 1:N
        [h(hI,wI),p(hI,wI)] = kstest2(X{hI},X{wI},'Alpha',alpha);
    end
end

end