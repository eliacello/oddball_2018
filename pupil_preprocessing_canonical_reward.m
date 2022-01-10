function [pupil_size, timestamps, sampling_factor, pSamplesRejected] = pupil_preprocessing_canonical_reward(timestamps, pupil, gx, gy, do_mad_outlier, do_fixationCorrection, do_bandpass_filter, do_interpolation, do_resampling, do_zscore)
% function [pupil_size, timestamps, sampling_factor, pSamplesRejected] = pupil_preprocessing_canonical_v03(timestamps, pupil_size, valid_recording, do_mad_outlier, do_bandpass_filter, do_interpolation, do_resampling, do_zscore)
%
%  Preprocessing of pupil recordings.
%  On a timeseries of pupil diameter, the following preprocessing steps can
%  be performed (call them with input to function):

%  1) Calculate the recorded sampling rate
%  2) 'do_mad_outlier':     Remove outliers outside 3 MAD
%  3) 'do_bandpass_filter': Bandpass filter
%  4) 'do_interpolation':   Interpolate bad states
%  5) 'do_resampling':      Resampling to 60Hz if needed
%  6) 'do_zscore':          Z-score
%
%  If you have recordings form both eyes, please preprocess each eye
%  separately, and average them afterwards.
%
%   Inputs
%   timestamps:      Vector of time stamps, in ms resolution, size: nx1
%   pupil_size:      Vector of pupil diameter, size: nx1
%   valid_recording: Vector indicating good and bad states of recording (signal lost,
%                    fixations outside your stimuli on the screen, blinks etc.)
%                    The results should be a vector, indicating for each time sample
%                    if it is valid (1) or not (0).  Samples marked with 0 will be interpolated
%                    in pupil_size later. Size: nx1
%
%   Outputs
%   pupil_size:       Preprocessed pupil_size vector.
%   timestamps:       Timestamps after interpolation, downsampling etc.
%   sampling_factor:  How much we did up- or downsample
%   pSamplesRejected: What proportion of samples got rejected?
%
%   Author: Antonius Wiehler <antonius.wiehler@gmail.com>
%   Original: 2018-03-16
%   Modified: 2018-10-25



%% PARAMETERS - PLEASE DO NOT TOUCH :)
% =========================================================================
blink           = 0;      % bad data in the pupil vector will be replaced by this
samplingrateFinal    = 50;     % sampling rate we want to have after preprocessing
blinkwindow     = 0.1;    % how many seconds before and after blink do we want to remove?
mad_cutoff      = 3;      % ouliers of how many mad should be rejected? (median deviance)
highpass_filter = 1/128;  % pupil highpass filter in Hz.
lowpass_filter  = 1;      % Filter everything faster than this Hz.
border_size=0.05;
ScreenX= 1400;
ScreenY= 1050; 

%% CHECKS
% =========================================================================
if nargin < 8
    error('Not enough input arguments.');
end

% correct vector orientation
timestamps = timestamps(:);
pupil_size = pupil(:);
% valid_recording = valid_recording(:);

if size(pupil_size) ~= size(timestamps)
    error('Inputs do not have the same length.');
end

%         & size(pupil_size) == size(valid_recording))

%% PREPROCESSING
% =========================================================================

% calculate sampling rate
% -------------------------------------------------------------------------
samplingrate = 500; %to change to 1000 for controls_S131 and controls_S181


% remove bad states
% -------------------------------------------------------------------------
% pupil_size(~valid_recording) = blink;


% exclude samples that are outside mad_cutoff - remove outliers
% -------------------------------------------------------------------------
if do_mad_outlier
    pupil_size = remove_outliers_mad_v01(pupil_size, mad_cutoff, blink);
end



% Delete pupil recording when fixation is at the perepheric of the screen.
% -------------------------------------------------------------------------
if do_fixationCorrection
    [pupil_size, outside_samples] = delete_periphical_fixations_v01(pupil_size, gx, gy, ScreenX, ScreenY, border_size, blink);
end 

% how many samples have been rejected?
% -------------------------------------------------------------------------
pSamplesRejected = sum(pupil_size == blink) ./ length(pupil_size);

% interpolate blinks and bad states +- window
% -------------------------------------------------------------------------
if do_interpolation
    [pupil_size, blink_indx] = interpolate_blinks_nonlinear_v01(pupil_size, blink, samplingrate, blinkwindow);
end

% band-pass filter to remove slow fluctuations and to smooth at the same
% time
% -------------------------------------------------------------------------
if do_bandpass_filter
    pupil_size = bandpass_filter_v02(pupil_size, samplingrate, highpass_filter, lowpass_filter);
end

% downsampling
% -------------------------------------------------------------------------
sampling_factor = round(samplingrate / samplingrateFinal);  % by what factor do we need to reduce?

if do_resampling
    if sampling_factor > 1.1
        
        % downsamples with low pass filtering, used for measures
        pupil_size = decimate(pupil_size, sampling_factor);
        
        % downsample without low pass filtering, used for markers etc.
        timestamps = downsample(timestamps, sampling_factor);
        
        fprintf('Time series was downsampled from %.2fHz to %.2fHz.\n', samplingrate, samplingrateFinal);
        
    elseif sampling_factor < 0.9
        
        x          = 1 : length(timestamps);
        xq         = 1 : sampling_factor : length(timestamps) + 1;
        timestamps = interp1(x, timestamps, xq)';
        pupil_size = interp1(x, pupil_size, xq)';
        timestamps = timestamps(1 : end - 2, :);  % remove last to samples because they are likely nan
        pupil_size = pupil_size(1 : end - 2, :);  % remove last to samples because they are likely nan
        
        fprintf('Time series was upsampled from %.2fHz to %.2fHz.\n', samplingrate, samplingrateFinal);
        
    else
        % if we are in the accuracy range, we do not touch
        sampling_factor = nan;  % just to return a value
    end
end

% z-score pupil time series
% -------------------------------------------------------------------------
if do_zscore
    pupil_size = zscore(pupil_size);
end

% give feedback on console
% -------------------------------------------------------------------------
fprintf('Pupil preprocessing done - %.2f%% of samples rejected.\n', sum(pSamplesRejected) * 100);  % done :)


end  % main function
