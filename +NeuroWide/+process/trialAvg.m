function avg = trialAvg(data, idx, tBef, tAft, varargin)

p = inputParser;
addParameter(p,'dim',1);
addParameter(p,'rem_baseline',0);

parse(p,varargin{:});

%%

dims = size(data);
dims(p.Results.dim) = tBef + tAft + 1;

avg = zeros(dims);
N = numel(idx);

for i = 1:N
    t = idx(i) - tBef : idx(i) + tAft;
    rem_t = idx(i) - tBef : idx(i) - 1;

    if p.Results.rem_baseline
        if p.Results.dim == 3
            rem = mean(data(:,:,rem_t), 3);
        else
            rem = mean(data(rem_t,:));
        end
    else
        rem = 0;
    end

    if p.Results.dim == 3
        avg = avg + 1 / N * (data(:,:,t) - rem);
    else
        avg = avg + 1 / N * (data(t,:,:) - rem);
    end
end

end