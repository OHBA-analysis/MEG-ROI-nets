function [Ttmp, ptmp, corrptmp, COPE] = perform_glm_with_randomise(edges, designMatrix, contrasts, standardise)
%PERFORM_GLM_WITH_RANDOMISE 
%
% [T, P, CORRP] = PERFORM_GLM_WITH_RANDOMISE(DATA, X, CONTRASTS, STANDARDISE)
%   runs a GLM of DATA = X Beta + noise. Returns T-statistics, uncorrected
%   and permutation-corrected p-values associated with the regression. 
%
% [T, P, CORRP, COPE] = ... also returns the contrast of parameter
%   estimates from the regression. 
%
% Randomise actually returns 1-p values in p and corrp.
%
% Note: randomise automatically demeans, so this is not helpful for finding
% mean effect over all subjects / sessions

nSessions = ROInets.cols(edges);

[checkMe, nEVs] = size(designMatrix);
assert(checkMe == nSessions, ...
       [mfilename ':BadDesign'], ...
       'Design matrix must have as many rows as the data have columns. \n');
   
[nContrasts, checkMe] = size(contrasts);
assert(checkMe == nEVs, ...
       [mfilename ':BadContrasts'], ...
       'Contrasts must have as many columns as EVs in the design matrix. \n');

if nargin < 4 || ~exist('standardise', 'var'),
    standardise = false;
else
    assert(islogical(standardise), ...
           [mfilename ':BadStandardise'], ...
           'Standardise input must be a logical value. \n');
end%if

%% Save out edges into nifti
resultsDir = tempdir;
inputNifti = fullfile(resultsDir, 'network_edges.nii.gz');

% check for NaNs - treat as subjects to be ignored. This is a good
% representation for missing data. 
%
% What will we do with NaNs? Maybe impute the value as the group mean?
% Maybe EM imputation?
% Instead, let's remove them outright. 
badSubjects = any(isnan(edges),1);
cleanEdges  = edges;
cleanEdges(:,badSubjects) = [];

for iS = ROInets.cols(cleanEdges):-1:1,
    formattedEdges(:,1,1,iS) = cleanEdges(:,iS);
end
save_avw(formattedEdges, inputNifti, 'f', [1 1 1 1]);
Ci = onCleanup(@() delete(inputNifti));


%% Construct design matrix
designMatrix(badSubjects,:) = [];
if standardise, 
    % demean and variance normalise
    X = bsxfun(@rdivide, bsxfun(@minus, designMatrix, mean(designMatrix)), ...
                         std(designMatrix));
else
    X = designMatrix;
end%if

% save out
designFile = fullfile(resultsDir, 'univariate_edge_test_design.mat');
ROInets.save_vest(X, designFile);
Cd = onCleanup(@() delete(designFile));

%% Construct contrasts
contrastFile = fullfile(resultsDir, 'univariate_edge_test_design.con');
ROInets.save_vest(contrasts, contrastFile);
Cc = onCleanup(@() delete(contrastFile));

%% Run randomise
outputNifti = fullfile(resultsDir, 'univariate_edge_test');

% call to randomise
command = sprintf('randomise -i %s -o %s -d %s -t %s -x --norcmask', ...
                  inputNifti, outputNifti, designFile, contrastFile);
              
% submit to terminal
ROInets.call_fsl_wrapper(command);

Co = onCleanup(@() delete([outputNifti '*.nii.gz']));

%% Produce nice COPEs for each edge 
% a cope is the difference in mean Z-converted correlations between each
% group
pinvxtx = pinv(designMatrix' * designMatrix);
pinvx   = pinvxtx * designMatrix';

for iEdge = ROInets.rows(edges):-1:1,
    COPE(iEdge,:) = contrasts * pinvx * edges(iEdge, :).';
end%for

%% Retrieve results
for iCon = nContrasts:-1:1,
   TstatFile{iCon}     = [outputNifti '_tstat' num2str(iCon) '.nii.gz'];
   pFile{iCon}         = [outputNifti '_vox_p_tstat' num2str(iCon) '.nii.gz'];
   corrpFile{iCon}     = [outputNifti '_vox_corrp_tstat' num2str(iCon) '.nii.gz'];
   
   Ttmp(:,iCon)     = read_avw (TstatFile{iCon});
   ptmp(:,iCon)     = read_avw (pFile{iCon});
   corrptmp(:,iCon) = read_avw (corrpFile{iCon});
end%for
end%perform_glm_with_randomise
% [EOF]