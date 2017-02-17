function surrogate = generate_phase_surrogates(data, preserve_correlation,use_box_cox)
%% Generate surrogate data though method 3 in Hurtado et al 2004. Statistical
% method for detection of phase locking episodes in neural oscillations. J
% Neurophysiol. 10.1152/jn.00853.2003.
%
% This scrambles the phase spectrum whilst preserving the amplitude spectrum
% The surrogate data is rescaled to the same mean and standard deviation
%
% surrogate = generate_phase_surrogate(data, preserve_correlation, mean_term, std_term)
% 
% INPUTS
%    - data is the original timeseries, [nsamples x nchannels]
%    - preserve_correlation - apply the same phase offset to each channel. False by default
%	 - use_box_cox. false by default. If this is true, then boxcox1.m from OSL will be used to 
%	   transform the data. Each channel gets its own transformation. The inverse transform will 
%	   automatically be applied. Offsets are applied prior to doing the transform and inverse
%	   transform, following Adam Baker's implementation. However, the final rescaling of the mean
% 	   could lead to the surrogate containing negative values, even if the original data does not.
%	   This should be fine for checking correlations, since none of the rescalings affect the correlation
%	   except that if a subsequent analysis assumes the surrogate envelopes are positive, this may fail
% 	   unless an additional offset is applied
%
% OUTPUTS
%    - surrogate - the surrogate time series, same size as data
% 
% Romesh Abeysuriya Nov 2016
% Andrew Quinn August 2015
% Adjusted by GC to improve speed of fft

if nargin < 3 || isempty(use_box_cox) 
	use_box_cox = false;
end

if nargin < 2 || isempty(preserve_correlation) 
    preserve_correlation = false; % Generate a different phase offset vector for each channel
end

[nSamples, nVars] = size(data);
mean_term = mean(data);
std_term = std(data,1); % Use N-1 standard deviation to match original code

if use_box_cox
	if ~exist('boxcox1')
		error('boxcox1.m not found. Has OSL been added to the path?');
	end

	lambda = nan(nVars,1);
	c = nan(nVars,1);
	for j = 1:nVars
		[~,lambda(j),c(j)] = boxcox1(data(:,j),0);
		data(:,j) = (data(:,1)+c(j)).^lambda(j);
	end
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

% If we rescaled the timeseries, undo the transformation, following 
% Adam's thesis where raising to the power of 1/lambda covers all cases
if use_box_cox
	for j = 1:nVars
		surrogate(:,j) = surrogate(:,j) - min(surrogate(:,j)) + min(data(:,j));
		surrogate(:,j) = surrogate(:,j).^(1/lambda(j));
	end
end

% demean
surrogate = bsxfun(@minus, surrogate, mean(surrogate));

% Normalise time_series
surrogate = bsxfun(@times, surrogate, std_term./std(surrogate,1)); 

% Reset the mean
surrogate = bsxfun(@plus, surrogate, mean_term);

end%generate_phase_surrogates