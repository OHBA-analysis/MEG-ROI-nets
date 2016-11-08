function surrogate = generate_phase_surrogate(data, preserve_correlation, mean_term, std_term)
%% Generate surrogate data though method 3 in Hurtado et al 2004. Statistical
% method for detection of phase locking episodes in neural oscillations. J
% Neurophysiol. 10.1152/jn.00853.2003.
%
% This scrambles the phase spectrum whilst preserving the amplitude spectrum
%
% surrogate = generate_phase_surrogate(data, preserve_correlation, mean_term, std_term)
% 
% INPUTS
%    - data is the original timeseries, [nsamples x nchannels]
%    - preserve_correlation - apply the same phase offset to each channel. False by default
%    - mean_term - Set the mean of the surrogate. By default, channels will have the same 
%      mean as the input data. Can be 1x1 or 1 x nchannels
%    - std_term - Set the standard deviation of each channel independently. Can be
%      1x1 or 1 x nchannels. By default, each channel is rescaled to the same standard deviation
%      as the input data channel
%
% OUTPUTS
%    - surrogate - the surrogate time series, same size as data
% 
% Romesh Abeysuriya Nov 2016
% Andrew Quinn August 2015
% Adjusted by GC to improve speed of fft

[nSamples, nVars] = size(data);

if nargin < 4 || isempty(std_term)
    % If the output std isn't defined, we'll use the std of the data
    std_term = std(data,1); % normalize by N, not N-1, to match original code
end

if nargin < 3 || isempty(mean_term)
    % Set surrogate mean to the same as the input data
    mean_term = mean(data);
end

if nargin < 2 || isempty(preserve_correlation) 
    preserve_correlation = false; % Generate a different phase offset vector for each channel
end

% get amplitude spectrum
amp_spec = fft(data, nSamples,1);      % no need to take abs as we can just add random phases

% Generate phase noise
% first element of spectrum is DC component, subsequent elements mirror
% each other about the Nyquist frequency, if present
n_components = floor((nSamples-1)/2); % Number of components that *aren't* DC or Nyquist
if preserve_correlation
    noise = bsxfun(@times,ones(1,nVars),rand(n_components,1).* (2 * pi)); % Same phase shift for each channel
else
    noise = rand(n_components, nVars) .* (2 * pi);
end

if rem(nSamples,2) % If odd number of samples, then Nyquist frequency is NOT present
    newPhase = [zeros(1,nVars); 1i .* noise; -1i .* flipud(noise)]; % second half uses conjugate phases, mirrored in order
else % Otherwise, include zero phase shift for the Nyquist frequency
    newPhase = [zeros(1,nVars); 1i .* noise; zeros(1,nVars); -1i .* flipud(noise)]; % second half uses conjugate phases, mirrored in order
end

% Make new phase scrabled spectrum
rand_spec = exp(newPhase).*amp_spec;

% Create new time_course
surrogate = ifft(rand_spec, nSamples);

% demean
surrogate = bsxfun(@minus, surrogate, mean(surrogate));

% Normalise time_series
surrogate = bsxfun(@times, surrogate, std_term./std(surrogate,1)); 

% Reset the mean
surrogate = bsxfun(@plus, surrogate, mean_term);

end%generate_phase_surrogates