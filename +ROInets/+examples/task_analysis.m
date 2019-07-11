function correlationMats = task_analysis(varargin)
%EXAMPLE  example task analysis for epoched data from several sessions
%
% correlationMats = EXAMPLE_TASK_ANALYIS()
% for help, type
% `help run_network_analysis'

%	Copyright 2015 OHBA
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


%	$LastChangedBy$
%	$Revision$
%	$LastChangedDate$
%	Contact: giles.colclough@gmail.com
%	Originally written on: MACI64 by Giles Colclough, 12-Jan-2015 15:25:47


% choose a binary ROI map. Take care that the resolution of the nifti file
% matches that of the source reconstruction.
parcellationDirectory = '/path/to/parcellation/';
parcelFile = fullfile(parcellationDirectory, ...
                      'parcellationFile.nii.gz');
                  
% choose a results directory
outDir = '/path/to/results/';

% set objects in and names for each session
dataFiles = {'/Users/gilesc/data/MND-malcolm/Motor_beta.oat/concatMfsession12_spm_meeg.mat'; ...
             '/Users/gilesc/data/MND-malcolm/Motor_beta.oat/concatMfsession120_spm_meeg.mat'; ...
             '/Users/gilesc/data/MND-malcolm/Motor_beta.oat/concatMfsession121_spm_meeg.mat'; ...
             '/Users/gilesc/data/MND-malcolm/Motor_beta.oat/concatMfsession122_spm_meeg.mat'};

         % setup the ROI network settings
Settings = struct();
Settings.spatialBasisSet          = parcelFile;                     % a binary file which holds the voxel allocation for each ROI - voxels x ROIs
Settings.gridStep                 = 8; % mm                         % resolution of source recon and nifti parcellation file
Settings.timeRange                = [0.01 3.99];                    % range of times to use for analysis
Settings.Regularize.do            = true;                           % use regularization on partial correlation matrices using the graphical lasso. 
Settings.Regularize.path          = 0.001;                          % This specifies a single, or vector, of possible rho-parameters controlling the strength of regularization. 
Settings.Regularize.method        = 'Friedman';                     % Regularization approach to take. {'Friedman' or 'Bayesian'}
Settings.Regularize.adaptivePath  = false;                          % adapt the regularization path if necessary
Settings.leakageCorrectionMethod  = 'closest';                      % choose from 'closest', 'symmetric', 'pairwise' or 'none'. 
Settings.nEmpiricalSamples        = 1;                              % convert correlations to standard normal z-statistics using a simulated empirical distribution. This controls how many times we simulate data of the same size as the analysis dataset
Settings.EnvelopeParams.windowLength = 1/40; % s                    % sliding window length for power envelope calculation. See Brookes 2011, 2012 and Luckhoo 2012. 
Settings.EnvelopeParams.useFilter = true;                           % use a more sophisticated filter than a sliding window average
Settings.EnvelopeParams.takeLogs  = true;
Settings.frequencyBands           = {[]};                           % a set of frequency bands for analysis. Set to empty to use broadband. The bandpass filtering is performed before orthogonalisation. 
Settings.timecourseCreationMethod = 'PCA';                          % 'PCA', 'mean', 'peakVoxel' or 'spatialBasis'
Settings.outputDirectory          = outDir;                         % Set a directory for the results output
Settings.groupStatisticsMethod    = 'mixed-effects';                % 'mixed-effects' or 'fixed-effects'
Settings.FDRalpha                 = 0.05;                           % false determination rate significance threshold
Settings.sessionName              = {'sess1', 'sess2', 'sess3', 'sess4'}; 
Settings.SaveCorrected.timeCourse    = false;
Settings.SaveCorrected.envelopes      = false;
Settings.SaveCorrected.variances     = false;
Settings.SaveCorrected.ROIweightings = false;
Settings.SubjectLevel.conditionLabel   = {'longvalidR', 'longvalidL', 'ShortValidRight', 'ShortValidLeft'};
Settings.SubjectLevel.designSummary    = {[1 0 0 0]', [0 1 0 0]', [0 0 1 0]', [0 0 0 1]'}; % summarise the design matrix by condition label
Settings.SubjectLevel.contrasts        = {[1 0 1 0]; [1 -1 1 -1]};  % each contrast is a new cell
Settings.GroupLevel.designMatrix       = [1 0                       % nSubjects x nEVs
                                          1 0
                                          0 1
                                          0 1];
Settings.GroupLevel.contrasts          = [1  1;  % contrast 1
                                          1 -1]; % contrast 2
% run the ROI network analysis
correlationMats = ROInets.run_network_analysis(dataFiles, Settings);
end%example_many_subj
% [EOF]
