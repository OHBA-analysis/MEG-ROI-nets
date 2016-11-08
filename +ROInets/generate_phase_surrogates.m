function surrogates = generate_phase_surrogates(data, nSurrogates, mean_term, std_term)
%% Generate surrogate data though method 3 in Hurtado et al 2004. Statistical
% method for detection of phase locking episodes in neural oscillations. J
% Neurophysiol. 10.1152/jn.00853.2003.
%
% This scrambles the phase spectrum whilst preserving the amplitude spectrum
%
% data is [nsamples x 1]
% n is an integer denoting the number of surrogates
%
% Andrew Quinn August 2015
% Adjusted by GC to improve speed of fft

[nSamples, nVars] = size(data);
assert(1==nVars,                          ...
      [mfilename ':IncorrectDataFormat'], ...
      'Data should be nSamples x 1. \n');

if nargin < 4 || isempty(std_term)
    % If the output std isn't defined, we'll use the std of the data
    % use inline std
    std_term = sqrt(sum(abs(bsxfun(@minus, data, sum(data)./nSamples)).^2)./nSamples);   
end

if nargin < 3 || isempty(mean_term)
    % Set surrogate mean to zero
    mean_term = 0;
end

% get amplitude spectrum
nfft     = nSamples + mod(nSamples,2); % speed improvements with nfft multiple of 2, or even a mutiple of 2. 
amp_spec = fft(data, nfft);      % no need to take abs as we can just add random phases

% Generate phase noise
% first element of spectrum is DC component, subsequent elements mirror
% each other about nfft/2 + 1
noise    = rand(nfft/2, nSurrogates) .* (2 * pi);
newPhase = [zeros(1,nSurrogates); ... % DC term has no phase
            1i .* noise;          ... % new phases, incl. central frequency
            -1i .* flipud(noise(1:end-1,:))]; % second half uses conjugate phases, mirrored in order

% Make new phase scrabled spectrum
rand_spec = bsxfun(@times, exp(newPhase), amp_spec);

% Create new time_course
surrogates = real(ifft(rand_spec, nfft));

% remove padded zeros
surrogates(nSamples+1:end,:) = [];

% demean
surrogates = bsxfun(@minus, surrogates, sum(surrogates,1)./nSamples);

% Normalise time_series
surrogates = bsxfun(@times, surrogates, ...
                    std_term ./ sqrt(sum(abs(surrogates).^2,1)./nSamples)); % inline standard using zero mean

if mean_term,
    surrogates = surrogates + mean_term;
end%if
end%generate_phase_surrogates