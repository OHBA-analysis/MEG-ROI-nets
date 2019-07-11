function out = check_inputs(Settings)
%CHECK_INPUTS	Checks properties of Settings structure
%
% Settings = ROInets.check_inputs(Settings) checks inputs and establishes 
%   default settings


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
%	$Revision: 374 $
%	$LastChangedDate: 2015-01-12 17:30:23 +0000 (Mon, 12 Jan 2015) $
%	Contact: giles.colclough@gmail.com
%	Originally written on: MACI64 by Giles Colclough, 17-Mar-2014 18:23:47

%% Define object
Inputs               = inputParser;
Inputs.CaseSensitive = false;
Inputs.FunctionName  = 'run_network_analysis';
Inputs.StructExpand  = true;  % If true, can pass parameter-value pairs in a struct
Inputs.KeepUnmatched = false; % If true, accept unexpected inputs


%% validation functions
numericValidFcn     = @(x) (isnumeric(x) && isscalar(x) && ...
                            ~isnan(x)    && ~isinf(x));

probabilityValidFcn = @(x) (numericValidFcn(x) && 0 <= x && x <= 1);

% check time range
validTimeRangeCheck = @(t) validateattributes(t, {'numeric'},         ...
                                              {'vector', 'numel', 2}, ...
                                              Inputs.FunctionName, 'timeRange');
increasingCheck     = @(t) assert(t(2) > t(1),                      ...
								  [mfilename ':TimeNotIncreasing'], ...
								  'timeRange must be an increasing two-component vector. \n');
timeRangeFcn        = @(x) isempty(x)   ||                     ...
                           iscell(x)    ||                     ...
                           (isvector(x) && (2 == length(x)) && ...
                           (x(2) > x(1)));
 
% check nifti files
isParcellationFcn   = @(x) (ischar(x)  && (exist(x, 'file')          ||   ...
                                           exist([x '.nii'], 'file') ||   ...
                                           exist([x '.nii.gz'], 'file'))) ...
                           || isnumeric(x) || islogical(x);
                                       
% check chosen orthog method
isProtocolFcn = @(c) ischar(validatestring(c,                                 ...
                                           {'none', 'symmetric',              ...
                                            'closest', 'pairwise'},           ...
                                           Inputs.FunctionName,               ...
                                           'leakageCorrectionMethod'));
                                        
% check chosen time-course method
isMethodFcn = @(c) ischar(validatestring(c,                            ...
                                        {'PCA', 'mean',                ...
                                         'spatialBasis', 'peakVoxel'}, ...
                                        Inputs.FunctionName,           ...
                                        'timecourseCreationMethod'));
                                    
% check group stats method
isGroupStatsMethodFcn = @(c) ischar(validatestring(c,                   ...
                                                   {'mixed-effects',    ...
                                                    'fixed-effects'},   ...
                                                   Inputs.FunctionName, ...
                                                   'groupStatisticsMethod'));
                                               
% check session names
isSessionNameFcn = @(c) ischar(c) || (iscell(c) && all(cellfun(@ischar, c)));

%% Inputs
Inputs.addParamValue('spatialBasisSet', [], isParcellationFcn);

Inputs.addParamValue('outputDirectory', [], @ischar);

Inputs.addParamValue('timeRange', [], timeRangeFcn);

RegularizeDefault = struct('do',        false, ...
                           'method',       [], ...
                           'path',         [], ...
                           'adaptivePath', [], ...
                           'Prior', struct('a', [], 'b', []));
Inputs.addParamValue('Regularize', RegularizeDefault, @isstruct);

protocolDefault = 'closest';
Inputs.addParamValue('leakageCorrectionMethod', protocolDefault, isProtocolFcn);

empiricalSamplesDefault = 0;
Inputs.addParamValue('nEmpiricalSamples', empiricalSamplesDefault, numericValidFcn);

ARmodelOrderDefault = 0;
Inputs.addParamValue('ARmodelOrder', ARmodelOrderDefault, numericValidFcn);

EnvParamsDefault = struct('windowLength', 2,    ...
                          'overlap',      0.6,  ...
						  'takeLogs',     true, ...
                          'useFilter',    false);
Inputs.addParamValue('EnvelopeParams', EnvParamsDefault, @isstruct);

FirstLevelDefault = [];
Inputs.addParamValue('FirstLevel', FirstLevelDefault, @(s) isempty(s) || isstruct(s));

SubjectLevelDefault = [];
Inputs.addParamValue('SubjectLevel', SubjectLevelDefault, @(s) isempty(s) || isstruct(s));

GroupLevelDefault = [];
Inputs.addParamValue('GroupLevel', GroupLevelDefault, @(s) isempty(s) || isstruct(s));

freqBandDefault = {[]};
Inputs.addParamValue('frequencyBands', freqBandDefault, @iscell);

tcsMethodDefault = 'PCA';
Inputs.addParamValue('timecourseCreationMethod', tcsMethodDefault, isMethodFcn);

groupStatsMethodDefault = 'mixed-effects';
Inputs.addParamValue('groupStatisticsMethod', groupStatsMethodDefault, isGroupStatsMethodFcn);

FDRalphaDefault = 0.05;
Inputs.addParamValue('FDRalpha', FDRalphaDefault, probabilityValidFcn);

SaveCorrectedDef = struct('timeCourses',   false, ...
                          'envelopes',     true,  ...
                          'variances',     true,  ...
                          'ROIweightings', false);
Inputs.addParamValue('SaveCorrected', SaveCorrectedDef, @isstruct);

Inputs.addParamValue('sessionName', 'session1', isSessionNameFcn);
Inputs.addParamValue('gridStep',    8,          numericValidFcn);

%% Parse inputs
Inputs.parse(Settings);

% add any further useful options
Settings                     = Inputs.Results;
Settings.dateRun             = datestr(now, 31);
Settings.nSessions           = max(double(ischar(Settings.sessionName)), ...
                                   length(Settings.sessionName));
Settings.nFreqBands          = length(Inputs.Results.frequencyBands);

% check any required inputs
if any(strcmpi('spatialBasisSet', Inputs.UsingDefaults)), 
    error([mfilename ':ParcelFileNotFound'], ...
          'You must provide a parcellation file or spatial basis set. \n');
end%if

% catch any settings not yet tailored for
if strcmpi('ica', Inputs.Results.spatialBasisSet),
    error('Performing the network analysis using the ICA results within the OIL structure is not yet possible. \n');
end%if

% check time range parameters
if ~isempty(Inputs.Results.timeRange),
    if iscell(Inputs.Results.timeRange),
        assert(isequal(length(Inputs.Results.timeRange), Settings.nSessions), ... 
               [Inputs.FunctionName ':BadTimeRangeCellLength'],             ...
               ['%s: timeRange input must be empty, a time range as a ',    ...
                'two-component vector, or a cell array of two-component ',  ...
                'vectors. The cell array must have the same length as the', ...
                'number of subjects to do. \n'],                            ...
               Inputs.FunctionName);
        cellfun(validTimeRangeCheck, Inputs.Results.timeRange);
		cellfun(increasingCheck, Inputs.Results.timeRange);
	else
		% not a cell - should be a vector. Will error if not. 
        validTimeRangeCheck(Inputs.Results.timeRange);
		increasingCheck(Inputs.Results.timeRange);		
    end%if
end%if

% check enveloping parameters
EnvelopeParams = Inputs.Results.EnvelopeParams;
numFieldNames = {'windowLength', 'overlap'};
for iff = 1:2,
    fieldName = numFieldNames{iff};
    if isfield(EnvelopeParams, fieldName),
        val = EnvelopeParams.(fieldName);
        assert(numericValidFcn(val) && val >= 0,            ...
               [mfilename ':IncorrectEnvField:' fieldName], ...
               ['EnvelopeParams.' fieldName ' must be a non-negative scalar. \n']);
    else
        EnvelopeParams.(fieldName) = EnvParamsDefault.(fieldName);
    end%if
end%for
if isfield(EnvelopeParams, 'useFilter'),
        assert(islogical(EnvelopeParams.useFilter),        ...
               [mfilename ':IncorrectEnvField:useFilter'], ...
               'EnvelopeParams.useFilter must be a logical value. \n');
else
    EnvelopeParams.useFilter = EnvParamsDefault.useFilter;
end%if
if isfield(EnvelopeParams, 'takeLogs'),
        assert(islogical(EnvelopeParams.takeLogs),        ...
               [mfilename ':IncorrectEnvField:takeLogs'], ...
               'EnvelopeParams.takeLogs must be a logical value. \n');
else
    EnvelopeParams.takeLogs = EnvParamsDefault.takeLogs;
end%if
EnvelopeParams.useHanningWindow = false; % not an option at the moment due to a bug in osl

% check save request params
ff = fieldnames(SaveCorrectedDef);
SaveCorrected = Inputs.Results.SaveCorrected;
for iff = 1:length(ff),
    fieldName = ff{iff};
    if isfield(SaveCorrected, fieldName),
        assert(islogical(SaveCorrected.(fieldName)),                  ...
               [mfilename ':IncorrectSaveCorrectedField:' fieldName], ...
               ['SaveCorrected.' fieldName ' must be a logical value. \n']);
    else
        SaveCorrected.(fieldName) = SaveCorrectedDef.(fieldName);
    end%if
end%for

Settings.EnvelopeParams      = EnvelopeParams;                             % incorporate defaults
Settings.SaveCorrected       = SaveCorrected;                              % incorporate defaults

% check regularization parameters
Settings.Regularize = check_regularization(Inputs.Results.Regularize, ...
                                           Inputs.Results.leakageCorrectionMethod);

% check first level options
if ~isempty(Settings.FirstLevel),
    check_first_level(Settings.FirstLevel);
	nSubjects             = check_subject_level(Settings.SubjectLevel, Settings.nSessions);
	Settings.nSubjects    = nSubjects;
    Settings.paradigm     = 'task';
else
	nSubjects             = Settings.nSessions;
    Settings.paradigm     = 'rest';
end%if


% and group level options
check_group_level(Settings.GroupLevel, nSubjects);

%% Output
out = Settings;
end%check_inputs


function Reg = check_regularization(Reg, leakageCorrectionMethod)
%CHECK_REGULARIZATION

assert(isfield(Reg, 'do') && islogical(Reg.do), ...
       [mfilename ':IncorrectRegularizationSpecification'], ...
       'The Regularize input must be a structure with a logical field ''do''. \n');
if Reg.do,
    if strcmpi(leakageCorrectionMethod, 'pairwise'),
        Reg.do = false;
        warning([mfilename ':noRegularizationForPairwise'],                    ...
                ['It is not possible to compute regularized partial ',         ...
                 'correlation using the pairwise leakage correction method, ', ...
                 'yet. \n   The regularization has been turned off. \n']);
    else
        allowedRegMethods = {'Bayesian', 'Friedman'};
        assert((isfield(Reg, 'method')                            ...
                && any(strcmp(Reg.method, allowedRegMethods))),   ...
               [mfilename ':IncorrectRegularizationMethod'],      ...
               ['Regularize.method must be specified as either ', ...
                '''Bayesian'' or ''Friedman''. \n']);
           
           if strcmp(Reg.method, 'Bayesian'),
               isPriorFormatted = @(x) numericValidFcn(x) && x >= 0;
               isPriorSpecified = @(Reg) isstruct(Reg.Prior)             ...
                                  && all(isfield(Reg.Prior, {'a', 'b'})) ...
                                  && isPriorFormatted(Reg.Prior.a)       ...
                                  && isPriorFormatted(Reg.Prior.b);
               if isfield(Reg, 'Prior') && ~isempty(Reg.Prior),
                   assert(isPriorSpecified(Reg),                            ...
                         [mfilename, 'IncorrectRegularizationPrior'], ...
                         'Check the specification of Regularize.Prior. \n');
               else % use default Kerman neutral Gamma prior (2011)
                   fprintf(['  %s: Bayesian regularization selected. ',       ...
                            'Using default Kerman neutral Gamma hyperprior ', ...
                            'on regularization parameter. \n'], mfilename);
                   Reg.Prior = struct('a', 1/3, 'b', 0);
               end%if
                      
           else % Friedman
               assert(isfield(Reg, 'path') && ~isempty(Reg.path) && isnumeric(Reg.path), ...
                      [mfilename ':IncorrectRegularizationPath'],                        ...
                      'Please input one or more regularization parameters in Regularize.path. \n');
               assert(isfield(Reg, 'adaptivePath') && islogical(Reg.adaptivePath), ...
                      [mfilename ':IncorrectAdaptive'],                            ...
                      'Pleast specify a logical parameter in Regularize.adaptivePath. \n');  
               assert(all(Reg.path >= 0) && all(diff(Reg.path) > 0),      ...
                      [mfilename ':IncorrectRegularizationPathValues'],   ... 
                      ['Regression parameters must be non-negative and ', ...
                       'increasing in value. \n']);
           end%if
    end%if
end%if
end%check_regularization

function check_first_level(FirstLevel)
%CHECK_FIRST_LEVEL

% are we running task analysis?
if isempty(FirstLevel),
    return
end%if

assert(isfield(FirstLevel, 'designSummary') && iscell(FirstLevel.designSummary), ...
       [mfilename ':noDesignSummary'],       ...
       'SubjectLevel.designSummary not specified properly. \n');
   
assert(isfield(FirstLevel, 'conditionLabel') && iscell(FirstLevel.conditionLabel), ...
       [mfilename ':noConditionLabel'],       ...
       'SubjectLevel.conditionLabel not specified properly. \n');

assert(isfield(FirstLevel, 'contrasts') && iscell(FirstLevel.contrasts), ...
       [mfilename ':noConrasts'],       ...
       'SubjectLevel.contrasts not specified properly. \n');
   
end%check_first_level

function nSubjects = check_subject_level(SubjectLevel, nSessions)
%CHECK_SUBJECT_LEVEL

% are we running task analysis?
if isempty(SubjectLevel),
    return
end%if

assert(isfield(SubjectLevel, 'subjectDesign') && ismatrix(SubjectLevel.subjectDesign), ...
       [mfilename ':noSubjectDesign'],       ...
       'SubjectLevel.subjectDesign not specified properly. \n');
   
[checkMe, nSubjects] = size(SubjectLevel.subjectDesign);
assert(checkMe == nSessions, ...
       [mfilename ':BadSubjectDesign'], ...
       'Subject level design matrix must have as many rows as sessions. \n');
   
end%check_subject_level

function check_group_level(GroupLevel, nSessions)
%CHECK_GROUP_LEVEL

% are we running a group level analysis?
if isempty(GroupLevel),
    return
end%if

assert(isfield(GroupLevel, 'designMatrix') && ismatrix(GroupLevel.designMatrix), ...
       [mfilename ':noDesignMatrix'],       ...
       'GroupLevel.designMatrix not specified properly. \n');
   
assert(isfield(GroupLevel, 'contrasts') && ismatrix(GroupLevel.contrasts), ...
       [mfilename 'noGroupContrasts'], ...
       'GroupLevel.contrasts not specified properly. \n');
   
[checkMe, nEVs] = size(GroupLevel.designMatrix);
assert(checkMe == nSessions, ...
       [mfilename ':BadGroupDesign'], ...
       'Group level design matrix must have as many rows as subjects. \n');
   
[~, checkMe] = size(GroupLevel.contrasts);
assert(checkMe == nEVs, ...
       [mfilename ':BadContrasts'], ...
       'Group Contrasts must have as many columns as EVs in the design matrix. \n');

end%check_group_level
% [EOF]
