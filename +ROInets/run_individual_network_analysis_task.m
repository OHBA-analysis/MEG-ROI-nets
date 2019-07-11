function CorrMats = run_individual_network_analysis_task(D,        ...
                                                         Settings, ...
                                                         resultsSaveName, ...
                                                         iSession)
%RUN_INDIVIDUAL_CORRELATION_ANALYSIS runs a single session network analysis
%
% CORRMATS = RUN_INDIVIDUAL_NETWORK_ANALYSIS(SPMMEEG, SETTINGS, SAVENAME)
%   performs a network correlation analysis between ROIs on a single subject. 
%
%   SPMMEEG is an SPM MEEG object.
%     data: nvoxels x time MEG data for analysis for a single subject. It
%           is assumed that the number of voxels matches the number of voxels in
%           the ROI basis set, which is passed in as a nifti file (see below). 
%
%
%   RESULTS_SAVE_NAME sets the location where the individual session
%   correlation matrices will be saved to disc. 
%
%
% CORRMATS = RUN_INDIVIDUAL_NETWORK_ANALYSIS(VOXELDATAFILE, SETTINGS, SAVENAME)
%   uses the name of a .mat file on disk, VOXELDATAILE. This .mat file
%   should contain three variables:
%     data: nvoxels x time MEG data for analysis for a single subject. It
%           is assumed that the number of voxels matches the number of voxels in
%           the ROI basis set. 
%
%     time: time-course for the MEG data
%     sampleRateInHz: also required. A constant sampling rate is assumed.
%
%   Relevant input parameters are held in SETTINGS:
%      spatialBasisSet          - the name of a binary nifti file which 
%                                 holds the voxel allocation for each ROI
%                                 (each ROI allocation is a volume), or a
%                                 spatial basis set from ICA (each spatial
%                                 map is a volume)
%                               - alternatively, pass in a voxels x ROIs
%                                 matrix. 
%
%      gridStep                 - grid step in mm of spatialBasisSet nifti
%                                 file. 
%
%      timeRange                - two-component vector of times to
%                                 sub-select a range for analysis. Leave as
%                                 empty ([]) to use the entire range
%                                 (default behaviour). Also can be passed
%                                 as a cell array of two-component vectors,
%                                 one for each subject. 
%      
%      Regularize               - Structure controlling the use of
%                                 regularization to find partial
%                                 correlation matrices, using the graphical
%                                 lasso. It has fields:
%            .do            : true or false: controls use of
%                             regularization. [false]
%            .method        : 'Bayesian' or 'Friedman': use Wang (2012)'s
%                             Bayesian graphical lasso or Friedman (2007)'s
%                             graphical lasso
%            .path          : This specifies a single, or vector of, 
%                             possible rho-parameters controlling the
%                             strength of regularization.
%                             If a vector is selected, 10-fold cross-validation
%                             is used to pick the optimal rho-parameter
%                             for each subject. The corrected Akaike
%                             Information Criterion is used to select
%                             the optimal rho. 
%                             [Friedman only]
%            .adaptivePath  : true or false: adapt path if the best
%                             regularization parameter is on the edge of
%                             the path. Then finesse the path three times. 
%                             [Friedman only] 
%            .Prior         : structure with fields controlling the shape
%                             of the gamma(x; a, 1/b) hyperprior on the 
%                             regularization parameter. 
%                             If this field is not set, the default Kerman 
%                             (2011) neutral hyperprior is used (a=1/3, b=0)
%                             [Bayesian only]
%                .a - shape
%                .b - 1/scale
%           for references, type `help ROInets.run_correlation_analysis'. 
%
%           If the Bayesian method is selected and the nEmpiricalSamples
%           options is used (see below), each empirical sample will use a
%           different regularization parameter at random from the posterior
%           distribution of regularization parameters used as part of the
%           MCMC scheme. nEmpiricalSamples should be relatively high to
%           adequately sample the space of 8000 values. 
%
%      leakageCorrectionMethod  - choose from 'symmetric', 'closest', 
%                                 'pairwise' or 'none'.
%                                 The symmetric method orthogonormalises 
%                                 the ROI time-courses in an
%                                 all-to-all symmetric fashion. 
%                                 The closest method starts with the
%                                 symmetric method, then lifts the
%                                 normality constraint to iterate to the
%                                 closest orthogonal matrix. 
%                                 The pairwise method extends voxel-wise
%                                 leakage correction methods to the parcel
%                                 time-courses. Correlations between pairs
%                                 of nodes are the average of correlations
%                                 found when one node is regressed from the
%                                 other, and then visa-versa.
%                                 'None' applies no spatial leakage
%                                 correction. 
%
%      SubjectLevel             - holds parameters for first-level analysis
%                                  .designSummary: cell array with one
%                                                  element for each
%                                                  regressor, holding a
%                                                  vector with one value
%                                                  for each condition
%                                  .conditionLabel: string condition labels
%                                                   to relate to D.trials
%                                  .contrasts: cell array of vectors, each
%                                              of length nRegressors
%
%      nEmpiricalSamples        - convert correlations to standard normal 
%                                 z-statistics using a simulated empirical 
%                                 distribution. This controls how many 
%                                 times we simulate null data of the same 
%                                 size as the analysis dataset. 
%                                 Set to zero to make normal assumptions, 
%                                 and decrease runtime. 
%                                 For datasets with many ROIs, only a few
%                                 samples are required - O(10). 
%
%      ARmodelOrder             - We tailor the empirical data to have the 
%                                 same temporal smoothness as the MEG data.
%                                 Set an order of 0 to use normality
%                                 assumptions. 
%                                 An order of 1 should be ok, but you could
%                                 use partial autocorrelation on some of 
%                                 the voxel time-courses to check this 
%                                 assumption. (Check out the matlab function
%                                 aryule, and google for partial
%                                 autocorrelation.)
%
%      EnvelopeParams           - Structure controlling how time-courses
%                                 are enveloped and down-sampled. Has
%                                 fields:
%
%          .windowLength :        sliding window length for power envelope 
%                                 calculation. 2 s is a good value. 
%                                 See Brookes 2011, 2012 and Luckhoo 2012. 
%
%          .overlap      :        Overlap on the sliding window. Try 0.6?  
%
%          .useFilter    :        Set true or false to use a more 
%                                 sophisticated downsamping operation than 
%                                 a moving average to find the power envelope. 
%                                 The passed frequency is 1.0/envelopeWindowLength.
%                                 The window overlap parameter is then
%                                 irrelevant. 
%          .takeLogs     :        Take logarithms of power envelopes if TRUE
%                                 (recommended), to improve normality
%
%      frequencyBands           - cell array of frequency bands for analysis. 
%                                 Set to empty to use broadband. 
%                                 E.g. {[4 8], [8 13], []}
%                                 The bandpass filtering is performed 
%                                 before orthogonalisation. 
%
%      timecourseCreationMethod - 'PCA', 'peakVoxel' or 'mean' for binary 
%                                 parcellations. 
%                                 Sets whether an ROI time-course is 
%                                 created using the mean over all voxels in 
%                                 the ROI, or by choosing the coefficients 
%                                 of the principal component accounting for 
%                                 the majority of the ROI's variance, or by
%                                 selecting the voxel with the greatest
%                                 variance.
%
%                                 If a spatial basis set (e.g. from ICA) is
%                                 passed in, 'spatialBasis' can be set,
%                                 which uses the PCA method, and
%                                 incorporating the whole-brain weightings
%                                 of each spatial map. 
%
%      groupStatisticsMethod    - Choose 'mixed-effects' or 'fixed-effects'
%                                 for group-level analysis. The former
%                                 performs a t-test on the means of the
%                                 first-level z-stats; the latter performs
%                                 a z-test. The t-test may have less power
%                                 but is also less sensitive to
%                                 inaccuracies in the scaling of
%                                 first-level z-scores - it is preferred,
%                                 but may need a large number of subjects
%                                 to produce interpretable results.
%
%      outputDirectory          - Set a directory for the results output
%
%      sessionName              - string used to construct filenames for
%                                 saving outputs. 
%
%      SaveCorrected            - Structure determining what is saved to
%                                 disc. Has fields:
%
%          .timeCourses     : orthgonalised ROI time-courses
%          .envelopes       : corrected ROI envelopes
%          .variances       : variances in each ROI
%          .ROIweightings   : weights used in each ROI used to construct
%                             single time-course
%
%
%   The ROI time-courses, corrected for source spread, and their envelopes, 
%   are saved in 
%     fullfile(oil.ROInetworks.outputDirectory, 'corrected-ROI-timecourses')
%
%   Various correlation matrices are output in CORRELATIONMATS. The
%   variable is a cell array, with one cell for each frequency band in the
%   analysis. Each cell holds a structure with fields:
%      
%       correlation                             : correlation between 
%                                                 corrected time-courses 
%                                                 (expected to be null) 
%                                                 nROIs x nROIs x nSessions
%
%       envCorrelation                          : correlation between
%                                                 corrected BLP envelopes
%                                                 nROIs x nROIs x nSessions
%
%       envPartialCorrelation                   : partial correlation
%                                                 between BLP envelopes
%                                                 nROIs x nROIs x nSessions
%
%       envPartialCorrelationRegularized        : L1-regularized partial 
%                                                 correlation between BLP 
%                                                 envelopes
%
%       env_z                                   : z-stats for envelope
%                                                 correlations
%                                                 nROIs x nROIs x nSessions
%
%       env_z_partial                           : z-stats for envelope
%                                                 partial correlations
%                                                 nROIs x nROIs x nSessions
%
%       env_z_partial_reg                       : z-stats for envelope
%                                                 partial correlations
%                                                 nROIs x nROIs x nSessions
%
%       Regularization                          : Chosen rho-parameter for
%                                                 L1 regularization in each
%                                                 subject
%
%       ARmodel                                 : AR coefficients estimated
%                                                 from MEG data for each
%                                                 session
%
%       H0sigma                                 : Width of empirical H0
%                                                 null correlation z-stats
%
%   These are also saved to disk at: 
%   fullfile(Settings.outputDirectory, 'ROInetworks_correlation_mats.mat')
%   
%   The analysis settings are also saved in the output directory. 
%
%   The methodology used in this pipeline is set out in 
%   Colclough, G. L., Brookes, M., Smith, S. M. and Woolrich, M. W., "A
%   symmetric multivariate leakage correction for MEG connectomes,"
%   NeuroImage [in revision]. 



%	Copyright 2014 OHBA
%	This program is free software: you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation, either version 3 of the License, or
%	(at your option) any later version.
%	
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%	
%	You should have received a copy of the GNU General Public License
%	along with this program.  If not, see <http://www.gnu.org/licenses/>.


%	$LastChangedBy: giles.colclough@gmail.com $
%	$Revision: 263 $
%	$LastChangedDate: 2014-10-23 11:30:39 +0100 (Thu, 23 Oct 2014) $
%	Contact: giles.colclough@gmail.com
%	Originally written on: MACI64 by Giles Colclough, 17-Mar-2014 23:39:39

fprintf('%s: starting analysis. \n', mfilename);


%% Parse input settings
% fprintf('%s: checking inputs. \n', mfilename);
% Settings = ROInets.check_inputs(Settings);

% make results directory
ROInets.make_directory(Settings.outputDirectory);

% save settings
outputDirectory = Settings.outputDirectory;
save(fullfile(outputDirectory, 'ROInetworks_settings.mat'), 'Settings');

% load parcellation
parcelFile = Settings.spatialBasisSet;
if ~isempty(parcelFile) && ischar(parcelFile), 
    spatialBasis = nii_quickread(parcelFile, Settings.gridStep);
elseif (islogical(parcelFile) || isnumeric(parcelFile)) && ismatrix(parcelFile), % expecting vox x parcels
    spatialBasis = parcelFile;
else
    error([mfilename ':ParcelReadFail'], ...
          'Failed to read parcel file %s. \n', ...
          parcelFile);
end%if
clear parcelFile

% is this function being run as part of a loop?
if nargin < 4 || ~exist('iSession', 'var') || isempty(iSession),
    iSession    = 1;
    sessionName = Settings.sessionName;
else
    sessionName = Settings.sessionName{iSession};
end%if
doRegularize    = Settings.Regularize.do;
tcsMethod       = Settings.timecourseCreationMethod;

switch lower(Settings.leakageCorrectionMethod),
    case 'pairwise',
        % this is a bit of a fudge at the moment. We perform the
        % orthogonalisation and correlation steps together
        protocol = 'none';
        
    otherwise
        protocol = Settings.leakageCorrectionMethod;
end%switch

%% Extracting timecourses
% load in data
% assume spm meeg object
% parse any string inputs automatically
D                  = spm_eeg_load(D);
fprintf('%s: loading data from file %s. \n', mfilename, D.fname);
   
%% Reduce ROIs to save RAM
assert(isequal(ROInets.rows(spatialBasis), D.nchannels), ...
       [mfilename ':ROIDataMisMatch'],                   ...
       'Number of voxels in data and ROI basis set do not match. \n');
   
allROImask = logical(ROInets.row_sum(abs(spatialBasis)));
% remove irrelevant voxels with no weight in any basis vector
spatialBasis(~allROImask, :) = [];
   
%% Parse time information
% we cannot have differing numbers of good and bad samples in each trial,
% as this is a nightmare for controlling correlation analyses. 
% Enforce the same number of timepoints in each trial.

badSamples           = any(D.badsamples(:,:,:));
trialsWithBadSamples = squeeze(any(badSamples(:,:,:),2));
if any(trialsWithBadSamples),
    % some trials have bad samples
    warning([mfilename ':BadTrialsFound'],                                  ...
            '%d trials found with bad samples. These will be excluded. \n', ...
            sum(trialsWithBadSamples));
end%if

% select time range from user-defined input. No need to account for bad
% samples now.
[time, timeInd, timeRange] = time_range(D.time, ...
                                        Settings.timeRange, iSession);

%% Identify trials to use and set up design matrix
[designMat, goodTrials] = set_up_first_level(D, Settings, ...
                                             find(trialsWithBadSamples));
 
%% Perform separate analyses in each frequency band
for iFreq = Settings.nFreqBands:-1:1,
    fprintf(' Correlation analysis for frequency band %d of %d. \n', ...
            iFreq, Settings.nFreqBands);
        
    % set up a string to identify the band
    if isempty(Settings.frequencyBands{iFreq}), 
        bandName = 'broadband';
    else
        bandName  = sprintf('%d-%dHz', ...
                            Settings.frequencyBands{iFreq}(1), ...
                            Settings.frequencyBands{iFreq}(2));
    end%if
    
    %% bandpass filter the data
    if ~isempty(Settings.frequencyBands{iFreq}),
        % need to turn off source space montage first
        m     = D.montage('getnumber');
        Dfilt = spm_eeg_filter(struct('D', D.montage('switch', 0), ...
                                      'band', 'bandpass', ...
                                      'freq', Settings.frequencyBands{iFreq}));
        % easy cleanup
        cleanD = onCleanup(@() delete(Dfilt));
        
        % move back to source space
        Dfilt = Dfilt.montage('switch', m);
        Dfilt.save;
    else
        cleanD = [];
        Dfilt  = D;
    end%if
    
    %% find node time-courses and envelopes
    if strcmpi(Settings.leakageCorrectionMethod, 'pairwise'),
        error([mfilename ':UnsupportedPairwiseMethod'], ...
              'Pairwise correction not supported for trial-wise data. \n');
    end%if
    fprintf([' Finding time courses with method %s;\n', ...
             ' Removing source leakage with orthogonalisation method %s.\n'], ...
            tcsMethod, Settings.leakageCorrectionMethod);
        
    for iT = 1:length(goodTrials),
        iTrial = goodTrials(iT);
        % find node timecourses for each ROI
        nodeDataUnCorr = ROInets.get_node_tcs(Dfilt(find(allROImask),      ...
                                                    find(timeInd),iTrial), ...
                                              real(spatialBasis),          ...
                                              tcsMethod);                  %#ok<FNDSB>
        % applying symmetric orthogonalisation if necessary
        nodeData = ROInets.demean(ROInets.remove_source_leakage(ROInets.demean(nodeDataUnCorr,2), protocol),2);
        
        % take power envelopes
        [nodeEnv(:,:,iT), time_ds] = ROInets.envelope_data(nodeData,        ...    
                                                           time,            ...
                                                           Settings.EnvelopeParams); %#ok<ASGLU>
    end%for loop over trials
    
    % save power envelopes
    if Settings.SaveCorrected.envelopes,
        saveDir = fullfile(Settings.outputDirectory, ...
                           'corrected-ROI-timecourses', filesep);
        ROInets.make_directory(saveDir);
        saveFile = fullfile(saveDir,                                      ...
                            sprintf('%s_%s_ROI_envelope_timecourses.mat', ...
                                    sessionName, bandName));
        save(saveFile, 'nodeEnv', 'time_ds');
    end%if
        
    %% Run correlation analysis 
    % calculate correlation matrices. 
    CorrMats{iFreq} = ROInets.run_correlation_analysis([],           ...
                                                       nodeEnv,      ...
                                                       Settings.Regularize);
    
    
    % Use an empirical null to enable conversion to z-stats
    transformSurrogates = ~Settings.EnvelopeParams.takeLogs;
    RegParams           = struct('do', Settings.Regularize.do, ...
                                 'rho', CorrMats{iFreq}.Regularization.mean);
    sigma = ROInets.find_permutation_H0_distribution_width(nodeEnv,                    ...
                                                           Settings.nEmpiricalSamples, ...
                                                           RegParams,                  ...
                                                           transformSurrogates);
          
    CorrMats{iFreq}.H0Sigma = sigma;
    
    % Store session name
    CorrMats{iFreq}.sessionName = sessionName;
    CorrMats{iFreq}.timeWindow  = timeRange;
    
    % free up some memory
    clear nodeData nodeEnv
    
    %% conversion of correlations to z-stats
    fprintf(' Converting correlations to normal z-stats\n');
    CorrMats{iFreq} = ROInets.convert_correlations_to_normal_variables(CorrMats{iFreq}, ...
                                                                       sigma,      ...
                                                                       doRegularize);
    
    %% Run first-level GLM
    CorrMats{iFreq}.firstLevel = run_first_level_glm(CorrMats{iFreq},    ...
                                                     designMat,          ...
                                                     Settings.SubjectLevel.contrasts, ...
													 D.fname);
    % clean up filtered object
    delete(cleanD);
end%loop over freq bands

%% save results to disc to provide backup of completed session analyses
save(resultsSaveName, 'CorrMats');




%%% END OF FUNCTION PROPER %%%

end%run_individual_correlation_analysis
%--------------------------------------------------------------------------






%--------------------------------------------------------------------------
function FirstLevel = run_first_level_glm(CorrMats, designMat, contrasts, fileName)
%RUN_FIRST_LEVEL_GLM
%


% input checking
[nTrials, nRegressors] = size(designMat);
nContrasts             = length(contrasts);
[~, nModes, checkMe]   = size(CorrMats.envCorrelation_z);
assert(checkMe == nTrials,         ...
      [mfilename ':LostTrials'],   ...
      'Number of trials must match columns of design matrix. \n');

assert(iscell(contrasts),               ...
       [mfilename ':NonCellContrasts'], ...
       'Contrasts must be a cell array. \n');
assert(all(cellfun(@length, contrasts) == nRegressors), ...
       [mfilename ':BadContrastFormat'],                ...
       'All contrasts must have the same length as number of regressors. \n');
   
% make sure contrasts are formatted as a cell array of column vectors
useContrasts = cell(1,nContrasts);
for iContrast = 1:nContrasts,
    useContrasts{iContrast} = contrasts{iContrast}(:);
end%for
   
% Precompute some helpful things
XtX       = designMat' * designMat;
[RXtX, p] = chol(XtX);
if ~p,
	invXtX = RXtX \ (RXtX' \ eye(nRegressors));
	pinvX  = RXtX \ (RXtX' \ designMat');
	hasBadEVs    = false;
	badContrasts = false(nContrasts, 1);
else
	% design matrix was rank deficient
	% is that because we have missing information for certain trial types?
	badEVs    = all(0 == designMat);
	hasBadEVs = any(badEVs);
	if hasBadEVs,
		warning([mfilename ':MissingTrials'],                   ...
			    '%s: file %s is missing trials for %d EVs. \n', ...
				mfilename, fileName, sum(badEVs));
		badContrasts = logical(cellfun(@(C) any(C(badEVs)), useContrasts));
		invXtX = pinv(XtX);
		pinvX  = invXtX * designMat';
	else
		error([mfilename ':RankDeficientDesign'],                     ...
			  ['%s: the design matrix is rank deficient. ',           ...
			   'Check that you''ve specified your EVs sensibly. \n'], ...
			  mfilename);
	end%if
end%if

% declare memory
[rho, prho, prhoReg] = deal(zeros(nModes, nModes, nContrasts));

% run GLM on each edge
for i = 1:nModes,
    for j = i+1:nModes,
        rho(i,j,:) = glm_fast_for_meg(squeeze(CorrMats.envCorrelation_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useContrasts, 0);
        prho(i,j,:) = glm_fast_for_meg(squeeze(CorrMats.envPartialCorrelation_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useContrasts, 0);
								  
	    % fill in uninformative values with NaN.
		if hasBadEVs,
			rho(i,j,badContrasts)  = NaN;
			prho(i,j,badContrasts) = NaN;
		end%if
        if isfield(CorrMats, 'envPartialCorrelationRegularized_z'),
        prhoReg(i,j,1:nContrasts) = glm_fast_for_meg(squeeze(CorrMats.envPartialCorrelationRegularized_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useContrasts, 0); 
        prhoReg(i,j,badContrasts) = NaN;
        else
            prhoReg(i,j) = 0;
        end%if
    end%for
end%for

% symmetrise and reformat
for iContrast = nContrasts:-1:1,
    FirstLevel(iContrast).cope.correlation                   = rho(:,:,iContrast) + rho(:,:,iContrast)';
    FirstLevel(iContrast).cope.partialCorrelation            = prho(:,:,iContrast) + prho(:,:,iContrast)';
    FirstLevel(iContrast).cope.partialCorrelationRegularized = prhoReg(:,:,iContrast) + prhoReg(:,:,iContrast)';
end%for

end%run_first_level_glm






%--------------------------------------------------------------------------
function [designMat, goodTrials, trialID, nConditions] = set_up_first_level(D, Settings, excludeTrials)
%SET_UP_FIRST_LEVEL creates the design matrix and trial identifiers

nConditions = length(Settings.SubjectLevel.conditionLabel);

% hold the relevant indices for trials in each condition
for iCondition = nConditions:-1:1,
    tI = D.indtrial(Settings.SubjectLevel.conditionLabel{iCondition}, ...
                    'GOOD');
    trialInds{iCondition} = ROInets.setdiff_pos_int(tI, excludeTrials);
    % check we've found something                               
    if isempty(trialInds{iCondition}), 
        warning([mfilename ':EmptyCondition'], ...
                'No good trials found in %s for condition %s. \n', ...
                D.fname, Settings.SubjectLevel.conditionLabel{iCondition});
    end%if
end%for

% extract a list of all good, relevant trials
goodTrials = sort([trialInds{:}]);

% generate a set of IDs linking trials to condition number
trialID = zeros(length(goodTrials), 1);
for iCondition = nConditions:-1:1,
    trialID(ismember(goodTrials,trialInds{iCondition})) = iCondition;
end%for

% check design matrix size
assert(all(cellfun(@length, Settings.SubjectLevel.designSummary) == nConditions), ...
       [mfilename ':DesignMatrixSizeFault'],                                      ...
       ['The design matrix summary must, in each cell, contain a vector of the ', ...
        'same length as the number of conditions. \n']);
% use an OSL function to generate the subject-specific design matrix
designMat = oat_setup_designmatrix(struct('Xsummary', {Settings.SubjectLevel.designSummary}, ...
                                          'trialtypes', trialID));
end%set_up_first_level
    
%--------------------------------------------------------------------------
function [t, tI, tR] = time_range(time, timeRange, iSession)
%TIME_RANGE selects time range for analysis of each session
% TIME is a vector of times
% TIMERANGE is either a cell array of two-component vectors, a single
% two-component vector, or a null vector

if isempty(timeRange),
    % use the whole time range
    t = time;
    tI = true(size(time));
    tR = [];
else
    % subselect time range
    if iscell(timeRange),
        tR = timeRange{iSession};
    else
        tR = timeRange;
    end%if
    validateattributes(tR, {'numeric'}, ...
                       {'vector', 'numel', 2, 'nondecreasing'}, ... % can have negative times in task data. 
                       'time_range', 'timeRange', 2);
                   
    tI = (time <= tR(2)) & (time >= tR(1));
    t  = time(tI);
end%if
end%time_range



%--------------------------------------------------------------------------
% [EOF]
