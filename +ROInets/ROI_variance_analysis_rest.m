function [SubjectLevel, GroupLevel] = ROI_variance_analysis_rest(varianceFiles, subjectDesign, groupDesign, contrasts)
%ROI_VARIANCE_ANALYSIS_REST performs power analyis on subjects in resting
%  state
%
% [SUBJECTLEVEL, GROUPLEVEL] = ROI_VARIANCE_ANALYSIS_REST(FILES, S_X, G_X,C)
%
% takes in cell array of filenames FILES pointing to variance files. 
% There should be nSessions files. Each is a .mat file holding an 
% nROIs x 1 array of variances. The subject design matrix, S_X, should be 
% nSessions x nSubjects. (Use the identity matrix if there is only one
% session per subject.)
%
% The group design, G_X, should be nSubjects x nEVs. 
% Finally, the contrasts should be a matrix of nContrasts x nEVs.
%
% The outputs are SUBJECTLEVEL, a structure with parameters beta which are
% the mean log variance in each ROI for each subject, and GROUPLEVEL, which
% holds the beta values from the group-level regression, together with
% statistics associated with each contrast. 
%
% [...] = ROI_VARIANCE_ANALYSIS_REST(MAT, ...) takes in an nROIs by
% nSessions matrix of ROI variances.
%
% Note that analysis is performed on the logarithm of variances: effect
% sizes relate to % changes in variance. 
%
% This code uses randomise at the top level: the p-values for any
% whole-group mean effect will be uninformative. 
%
% Unlike in randomise's output, the p-values are p-values: a p-value of
% 0.001 is significant. (Not 1-p).

% parse input and transform to log-space
[variances, nSessions, nROIs] = parse_input(varianceFiles);
assert(all(variances(:) > 0), ...
	   'Expecting variances to be greater than 0. What''s up with your data?\n');
   
logVariances = log(variances);

% check dimensions
[checkMe, nSubjects] = size(subjectDesign);
assert(checkMe == nSessions, ...
	   'Subject design should be nSessions x nSubjects. \n');
[checkMe, nEVs] = size(groupDesign);
assert(checkMe == nSubjects, ...
	   'Group design should be nSubjects x nEVs. \n');
[nContrasts, checkMe] = size(contrasts);
assert(checkMe == nEVs, ...
	   'Contrasts should be nContrasts x nEVs. \n');

% perform subject-level glm
SubjectLevel              = struct();
SubjectLevel.design       = subjectDesign;
SubjectLevel.nSessions    = nSessions;
SubjectLevel.nSubjects    = nSubjects;
SubjectLevel.nROIs        = nROIs;
SubjectLevel.analysisTime = datestr(now);

SubjectLevel.beta         = subjectDesign \ logVariances';
SubjectLevel.meanVariance = exp(SubjectLevel.beta);

% perform group-level glm using randomise.
GroupLevel              = struct();
GroupLevel.design       = groupDesign;
GroupLevel.nSubjects    = nSubjects;
GroupLevel.nEVs         = nEVs;
GroupLevel.nContrasts   = nContrasts;
GroupLevel.analysisTime = datestr(now);

GroupLevel.beta               = groupDesign \ SubjectLevel.beta;
GroupLevel.parameterEstimates = contrasts * GroupLevel.beta;
[T, OneMinusp, corrOneMinusp] = ROInets.perform_glm_with_randomise(SubjectLevel.beta', groupDesign, contrasts, false);
for iContrast = nContrasts:-1:1,
	GroupLevel.contrast(iContrast).contrast = contrasts(iContrast,:);
	GroupLevel.contrast(iContrast).T             = T(:,iContrast);
	GroupLevel.contrast(iContrast).p             = 1 - OneMinusp(:,iContrast);
	GroupLevel.contrast(iContrast).minusLogp     = -log10(1 - OneMinusp(:,iContrast));
	GroupLevel.contrast(iContrast).corrp         = 1 - corrOneMinusp(:,iContrast);
	GroupLevel.contrast(iContrast).minusLogCorrp = -log10(1 - corrOneMinusp(:,iContrast));
end%for
end%ROI_variance_analysis_rest

function [V, nSessions, nROIs] = parse_input(varianceFiles)
%PARSE_INPUT extract variances from files

if iscell(varianceFiles),
	nSessions = length(varianceFiles);
	for iSession = nSessions:-1:1,
		tmp = load(varianceFiles{iSession});
		
		% this assumes that each file will have the same number of elements
		% stored
		V(:,iSession) = tmp.ROIvariances(:);
	end%for
	nROIs = ROInets.rows(V);
	
elseif ismatrix(varianceFiles),
	[nROIs, nSessions] = size(varianceFiles);
	V = varianceFiles;
	
else
	error([mfilename ':BadInput'], ...
		  'Expecting first input to be file list or a matrix of variances. \n');
end%if
end%parse_input
% [EOF]