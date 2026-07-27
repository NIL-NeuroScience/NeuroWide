function f = waterfall(data)

steps = prctile(data, [0.1, 99.9]);

steps = steps(2,2:end) - steps(1,1:end-1);
steps = cumsum(steps);
data = data - [0, steps];

f = figure;
plot(data);

end