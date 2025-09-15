function [p] = f_SRtest(X,alpha)

N = numel(X);

% h = zeros(N);
p = zeros(N);

for hI = 1:N
    for wI = 1:N
        [p(hI,wI)] = signrank(X{hI},X{wI},'alpha',alpha);
    end
end

end