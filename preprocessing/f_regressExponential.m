function video_out = f_regressExponential(video, brain_mask)
    % Function to regress an exponential fit from a 3D video using a brain mask
    %
    % Inputs:
    % video - 3D video (height x width x frames)
    % brain_mask - 2D mask (height x width) with NaN for non-brain regions and 1 for brain regions
    %
    % Output:
    % video_out - 3D video with the exponential fit regressed out

    % Ensure the brain mask is logical
    brain_mask = ~isnan(brain_mask);

    % Reshape the video into 2D (pixels x frames)
    [height, width, num_frames] = size(video);
    video_reshaped = reshape(video, [], num_frames);

    % Extract the brain region signals
    brain_signals = video_reshaped(brain_mask(:), :);

    % Calculate the average signal within the brain mask
    avg_signal = mean(brain_signals, 1, 'omitnan');

    % Define the exponential function to fit
    exp_func = @(p, t) p(1) * exp(-p(2) * t) + p(3);

    % Time vector
    t = (0:num_frames-1)';

    % Initial guess for the parameters [amplitude, decay rate, offset]
    p0 = [max(avg_signal) - min(avg_signal), 0.01, min(avg_signal)];

    % Fit the exponential function to the average signal
    opts = optimset('Display', 'off');
    p_fit = lsqcurvefit(exp_func, p0, t, avg_signal, [], [], opts);

    % Generate the fitted exponential signal
    fitted_signal = exp_func(p_fit, t);

    % Reshape the fitted exponential signal for regression
    % fitted_signal_reshaped = reshape(fitted_signal, 1, 1, []);

    % Initialize the output video
    video_out = zeros(size(video));

    % Perform linear regression for each pixel
    for i = 1:height
        for j = 1:width
            % Extract the pixel time series
            pixel_signal = squeeze(video(i, j, :));

            % Fit a linear model: pixel_signal ~ fitted_signal
            X = [fitted_signal, ones(num_frames, 1)]; % Design matrix
            beta = X \ pixel_signal; % Linear regression coefficients

            % Subtract the fitted exponential component
            video_out(i, j, :) = pixel_signal - X(:, 1) * beta(1);
        end
    end
end