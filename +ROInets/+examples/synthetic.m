function correlationMats = synthetic(outDir)
	% Example analysis using synthetic data (does not require an input meeg object
	% 
	% INPUTS
	% - outDir: Directory in which to write synthetic parcel timecourses and analysis output
	if nargin < 1 || isempty(outDir) 
		outDir = fullfile(pwd,'ROInets_synthetic_example');
		fprintf('Writing output to %s\n',outDir);
		mkdir(outDir)
	end
	
	% Write synthetic parcel timecourses
	dataFile = fullfile(outDir,'synthetic_data.mat');
	% Use an autoregressive model to generate smooth data
	% follow http://www.mathworks.co.uk/help/signal/examples/linear-prediction-and-autoregressive-modeling.html
	Fs       = 100; %Hz
	duration = 60; %s
	time     = 0:1.0/Fs:duration;
	nSamples = length(time);
	b        = fir1(1024, 0.5);
	nVoxels  = 3;
	[ARfilterTerms, ARnoiseVar] = lpc(b, 7); 
	% Generate data from a covariance matrix and smooth
	C    = [1  -0.1 0.6
	       -0.1   1 0.3
	        0.6 0.3   1] * ARnoiseVar;
	u    = chol(C)' * randn(nVoxels, nSamples);
	data = filter(1, ARfilterTerms, u.').';
	figure('Name', 'Input data', 'Color', 'w');
	plot(time.', data.');
	% Save to file
	sampleRateInHz = Fs;
	save(dataFile, 'data', 'time', 'sampleRateInHz');
	
	% Choose a binary ROI map. 
	spatialBasis = eye(3);
                 
	% set a save file name
	resultsName = fullfile(outDir, 'syntheticResults');
	sessionName = 'myBestSubject';

	% setup the ROI network settings
	Settings = struct();
	Settings.spatialBasisSet          = spatialBasis;                   % a binary file which holds the voxel allocation for each ROI - voxels x ROIs
	Settings.gridStep                 = 8; % mm                         % resolution of source recon and nifti parcellation file
	Settings.Regularize.do            = true;                           % use regularization on partial correlation matrices using the graphical lasso. 
	Settings.Regularize.path          = logspace(-9,2,80);              % This specifies a single, or vector, of possible rho-parameters controlling the strength of regularization. 
	Settings.Regularize.method        = 'Friedman';                     % Regularization approach to take. {'Friedman' or 'Bayesian'}
	Settings.Regularize.adaptivePath  = true;                           % adapth the regularization path if necessary
	Settings.leakageCorrectionMethod  = 'closest';                      % choose from 'closest', 'symmetric', 'pairwise' or 'none'. 
	Settings.nEmpiricalSamples        = 8;                              % convert correlations to standard normal z-statistics using a simulated empirical distribution. This controls how many times we simulate data of the same size as the analysis dataset
	Settings.ARmodelOrder             = 1;                              % We tailor the empirical data to have the same temporal smoothness as the MEG data. An order of 1 should be ok.
	Settings.EnvelopeParams.windowLength = 2; % s                       % sliding window length for power envelope calculation. See Brookes 2011, 2012 and Luckhoo 2012. 
	Settings.EnvelopeParams.useFilter    = true;                        % use a more sophisticated filter than a sliding window average
	Settings.EnvelopeParams.takeLogs  = true;
	Settings.frequencyBands           = {[13 30]};                      % a set of frequency bands for analysis. Set to empty to use broadband. The bandpass filtering is performed before orthogonalisation. 
	Settings.timecourseCreationMethod = 'spatialBasis';                 % 'PCA', 'peakVoxel' or 'spatialBasis'
	Settings.outputDirectory          = outDir;                         % Set a directory for the results output
	Settings.groupStatisticsMethod    = 'mixed-effects';                % 'mixed-effects' or 'fixed-effects'
	Settings.FDRalpha                 = 0.05;                           % false determination rate significance threshold
	Settings.sessionName              = sessionName; 
	Settings.SaveCorrected            = struct('timeCourses',   false, ...  % save corrected timecourses
	                                           'envelopes',     true,  ...  % save corrected power envelopes
	                                           'variances',     false);     % save mean power in each ROI before correction

	% run the ROI network analysis
	Settings        = ROInets.check_inputs(Settings);
	correlationMats = ROInets.run_individual_network_analysis(dataFile, Settings, resultsName);

	% Want to run an analysis on many subjects? Have a look at 
	% run_network_analysis to see the suggested steps. 

	% show results
	figure('Name', 'node correlation matrix', 'Color', 'w');
	imagesc(correlationMats{1}.envCorrelation);
	axis square
	colorbar

	% plot envelopes
	load(fullfile(outDir, 'corrected-ROI-timecourses', ...
	              'myBestSubject_13-30Hz_ROI_envelope_timecourses.mat'));
	figure('Name', 'envelope timecourses', 'Color', 'w');
	plot(time_ds', nodeEnv')
