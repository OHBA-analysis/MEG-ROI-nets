function correlationMats = run_network_analysis(Dlist, Settings)
%RUN_NETWORK_ANALYSIS	ROI correlation and partial correlation matrices
%
% [CORRELATIONMATS] = run_network_analysis(DLIST, Settings) runs a
%   correlation-based network analysis between ROIs on a set of
%   source-reconstructed data. These are provided as a cell array of SPM12
%   MEEG objects, or their filenames, in DLIST. 
%
%   Relevant input parameters are held in Settings:
%      spatialBasisSet          - the name of a binary nifti file which 
%                                 holds the voxel allocation for each ROI
%                                 (each ROI allocation is a volume), or a
%                                 spatial basis set from ICA (each spatial
%                                 map is a volume)
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
%      subjectsToDo             - vector of indices indicating which of the 
%                                 source reconstruction sessions to analyse
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
%      leakageCorrectionMethod  - choose from 'symmetric', 'pairwise' or 
%                                 'none'.
%                                 The symmetric method is recommended. It
%                                 orthogonalises the ROI time-courses in an
%                                 all-to-all symmetric fashion. 
%                                 The pairwise method extends voxel-wise
%                                 leakage correction methods to the parcel
%                                 time-courses. Correlations between pairs
%                                 of nodes are the average of correlations
%                                 found when one node is regressed from the
%                                 other, and then visa-versa.
%                                 'None' applies no spatial leakage
%                                 correction. 
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
%      timecourseCreationMethod - 'PCA' or 'peakVoxel'  for binary 
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
%  That's all you need for a resting-state analyis. 
%  If you want to run a task analysis on epoched data, you'll also need
%       FirstLevel             - holds parameters for first-level analysis
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
%       SubjectLevel           - holds parameters to deal with multiple
%                                sessions
%									.subjectDesign: nSessions x nSubjects
%									                matrix
%
%
%   In either case, if you want to run a group-level GLM analysis on each
%   network edge, correcting for family-wise error, you'll also need
%        GroupLevel              - holds parameter for group-level GLM
%                                  .designMatrix: nSessions x nRegressors
%                                                 matrix
%                                  .contrasts:    nContrasts x nRegressors
%                                                 matrix
%   There is a (slightly) different syntax at first, subject and group level. 
%
%
%   The ROI time-courses, corrected for source spread, and their envelopes, 
%   are saved in 
%     fullfile(Settings.outputDirectory, 'corrected-ROI-timecourses')
%
%   Various correlation matrices are output in CORRELATIONMATS. The
%   variable is a cell array, with one cell for each frequency band in the
%   analysis. Each cell holds a structure with fields:
%
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
% For a resting analysis, you'll also get:
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
%       envCorrelation_z                        : z-stats for envelope
%                                                 correlations
%                                                 nROIs x nROIs x nSessions
%
%       envPartialCorrelation_z                 : z-stats for envelope
%                                                 partial correlations
%                                                 nROIs x nROIs x nSessions
%
%       envPartialCorrelationRegularized_z      : z-stats for envelope
%                                                 partial correlaitons
%                                                 nROIs x nROIs x nSessions
%
%       groupEnvCorrelation_z                   : group result of z-test on
%                                                 mean of first-level
%                                                 z-statistics
%
%       groupEnvPartialCorrelation_z            : group result of z-test on
%                                                 mean of first-level
%                                                 z-statistics
%
%       groupEnvPartialCorrelationRegularized_z : group result of z-test on
%                                                 mean of first-level
%                                                 z-statistics
%   For a task analysis, you'll also get 
%       firstLevel: a struct of size equal to number of first level
%                   contrasts, with fields
%                       .cope : parameter estimates for each subject, held
%                               as stacked netmats with subjects in the
%                               third dimension
%                       .contrast : the contrast used
%                       .conditionLabel : the labels to which the contrast
%                                         refers
%
%   If you run a group GLM, you'll also get
%       groupLevel: a struct of size equal to number of first level
%                   contrasts, with fields
%                      .correlation
%                      .partialCorrelation
%                      .partialCorrelationRegularized
%                   Each contains univariate T-stats, the uncorrected
%                   p-values, the family-wise error p-values, and a
%                   significance matrix H under FDR correction. These
%                   results are held as netmats, with the group contrasts
%                   in the third dimension
%  
%       
%
%   These are results also saved to disk at: 
%   fullfile(Settings.outputDirectory, 'ROInetworks_correlation_mats.mat')
%   
%   The analysis settings are also saved in the output directory. 
%
%   References:
%   Colclough, G.L. et al., "A symmetric multivariate leakage correction
%   for MEG connectomes," NeuroImage 117, pp. 439-448 (2015). 



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


% JH, Jul-2019: removed this from the doc, because it does not seem to be implemented.
%
%      nParallelWorkers         - analyse several sessions in parallel, 
%                                 if non-zero. 
%                                 Specifies number of matlab workers to use. 
%                                 Caution: these tasks are very memory-heavy! 
%                                 You might slow down because you max out 
%                                 on RAM. 

%	$LastChangedBy$
%	$Revision$
%	$LastChangedDate$
%	Contact: giles.colclough@gmail.com
%	Originally written on: MACI64 by Giles Colclough, 17-Mar-2014 17:53:55

fprintf('%s: starting analysis. \n', mfilename);

%% Parse input settings
fprintf('%s: checking inputs. \n', mfilename);
Settings = ROInets.check_inputs(Settings);

assert(numel(Dlist) == Settings.nSessions, ...
       [mfilename ':WrongNoSessions'],     ...
       'Number of sessions should match number of D objects passed in. \n');

% make results directory
ROInets.make_directory(Settings.outputDirectory);

% save settings
outputDirectory = Settings.outputDirectory;
save(fullfile(outputDirectory, 'ROInetworks_settings.mat'), 'Settings');

%% Run correlation analysis on each subject

fprintf('%s: Running correlation analysis. \n', mfilename);

for iSession = Settings.nSessions:-1:1,
    fprintf('\n\n%s: Individual correlation analysis for file %d out of %d\n', ...
            mfilename, Settings.nSessions - iSession + 1, Settings.nSessions);

    D                          = Dlist{iSession};
    sessionName                = Settings.sessionName{iSession};
    matsSaveFileName{iSession} = fullfile(outputDirectory,                                      ...
                                          sprintf('%s_single_session_correlation_mats_tmp.mat', ...
                                                  sessionName));

    if strcmpi(Settings.paradigm, 'task'),
        mats{iSession} = ROInets.run_individual_network_analysis_task(D,                          ...
                                                                 Settings,                   ...
                                                                 matsSaveFileName{iSession}, ...
                                                                 iSession);
    elseif strcmpi(Settings.paradigm, 'rest'),
        mats{iSession} = ROInets.run_individual_network_analysis(D,                          ...
                                                                 Settings,                   ...
                                                                 matsSaveFileName{iSession}, ...
                                                                 iSession);
    else
        error([mfilename ':BadParadigm'], ...
              'Unrecognised paradigm %s. \n', Settings.paradigm);
    end%if
end%for

% reformat results - correlationMats is a cell array of frequency bands
correlationMats = ROInets.reformat_results(mats, Settings);

% save current results: will write over later
% just in case of crash at group level
saveFileName = fullfile(outputDirectory, 'ROInetworks_correlation_mats.mat');
save(saveFileName, 'correlationMats');
clear mats

%% Subject-level analysis to average over sessions in a fixed-effects manner
correlationMats = ROInets.do_subject_level_glm(correlationMats, Settings);

%% Group-level analysis
% Find whole group means
if strcmpi(Settings.paradigm, 'rest'),
    correlationMats = ROInets.do_group_level_statistics(correlationMats, Settings);
end%if

% Perform group-level GLM
if ~isempty(Settings.GroupLevel),
    correlationMats = ROInets.do_group_level_glm(correlationMats, Settings);
end%if

%% save matrices
fprintf('\n%s: Saving Results. \n', mfilename);

% save collected results
save(saveFileName, 'correlationMats');

% we stored individual results as we went along, in case of crash. Delete
% them if we've safely got to this stage. 
for iSession = 1:length(matsSaveFileName),
    delete(matsSaveFileName{iSession});
end%for
 
% tidy output of funciton
Settings.correlationMatsFile = saveFileName;
save(fullfile(outputDirectory, 'ROInetworks_settings.mat'), 'Settings');

fprintf('%s: Analysis complete. \n\n\n', mfilename);
end%run_network_analysis
% [EOF]
