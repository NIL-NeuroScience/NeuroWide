function [peak_inds] = f_peakFinder(roi, amp_cutoff, width, ploton)
%% finds peaks in timeseries
%excludes for motion
%excludes for amplitude

%inputs
    %roi = timeseries you'd like to extract peaks from (usually filtered)
    %t = corresponding time vector for roi
    %fd = corresponding motion vector for roi
    %amp_cutoff = instantaneous amplitude cutoff for peak detection
    %rng_bnds = range around peak  (vector from -range:range in samples)
        %example: rng_bnds=-20:20 would exclude peaks with motion in that
        %range around peaks
    %ploton = 1 if you want plots, 0 otherwise
    
%outputs
    %peak_times = time of peaks in s
    %peak_inds = index of peak in vector

%% peak identification

t = (1:numel(roi))/10;
fd = zeros(size(roi));

peak_times=[];
peak_inds=[];
tr=t(2)-t(1);

roi_peaks=movmax(roi, width);
imvmax = find(roi==roi_peaks);
mvmax = roi(imvmax);
time_peak=t(imvmax);

%% amplitude calc

hb = hilbert(roi);
amp = abs(hb);

%% loop through each peak, exclude for amplitude and motion
num_peaks=length(time_peak);

for i=1:num_peaks
    peak=imvmax(i);
    peak_t=time_peak(i);

    rng_bnds = 0;
    peak_rng=peak+rng_bnds;
    if (peak>abs(rng_bnds(1))) & (peak<length(t)-abs(rng_bnds(1)))
        if max(fd(peak_rng))<0.3
            if max(amp(peak+rng_bnds))>amp_cutoff
                peak_times=[peak_times peak_t];
                peak_inds=[peak_inds peak];
            end
        end
    end
end

if ploton==1
 figure(); plot(t, roi); hold on; plot(t(peak_inds), roi(peak_inds), '*')
end

end