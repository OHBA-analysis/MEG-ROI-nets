function correlationMats = many_subjects(Dlist,parcelFile,outDir,sessionName)
    % Example analysis for many subjects
    % 
    % INPUTS
    % - Dlist: array of meeg objects or filenames for meeg objects
    % - parcelFile: choose a binary ROI map. Take care that the resolution of the 
    %               nifti file matches that of the source reconstruction.
    % - outdir: choose a results directory
    % - sessionName: Optionally specify name of sessions (one for each entry in Dlist)

    if nargin < 4 || isempty(sessionName) 
        sessionName = arrayfun(@(x) sprintf('session_%d',x),1:length(Dlist),'UniformOutput',false);
    end

    % Set up the ROI network settings
    Settings = struct();
    Settings.spatialBasisSet          = parcelFile;                     % a binary file which holds the voxel allocation for each ROI - voxels x ROIs
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
    Settings.EnvelopeParams.takeLogs  = true;                           % perform analysis on logarithm of envelope. This improves normality assumption
    Settings.frequencyBands           = {[8 13], [13 30], []};          % a set of frequency bands for analysis. Set to empty to use broadband. The bandpass filtering is performed before orthogonalisation. 
    Settings.timecourseCreationMethod = 'spatialBasis';                 % 'PCA',  'peakVoxel' or 'spatialBasis'
    Settings.outputDirectory          = outDir;                         % Set a directory for the results output
    Settings.groupStatisticsMethod    = 'fixed-effects';                % 'mixed-effects' or 'fixed-effects'
    Settings.FDRalpha                 = 0.05;                           % false determination rate significance threshold
    Settings.sessionName              = sessionName; 
    Settings.SaveCorrected            = struct('timeCourses',   false, ...  % save corrected timecourses
                                               'envelopes',     true,  ...  % save corrected power envelopes
                                               'variances',     false);     % save mean power in each ROI before correction

    % Run the ROI network analysis
    correlationMats = ROInets.run_network_analysis(Dlist, Settings);
