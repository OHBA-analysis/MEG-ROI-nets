function [p, beta] = mean_connectivity_statistic_group_level(netmat, Settings, iContrast)
%MEAN_CONNECTIVITY_STATISTIC_GROUP_LEVEL
%
% [p, beta] = mean_connectivity_statistic_group_level(netmat, Settings, iContrast)
%
% Tests for a difference in mean connectivity over all edges
%
%  The model is set up by choosing a contrast from the group level
%  statistics contrasts in Settings. 
%
%  Netmat needs to contain an nNodes x nNodes x nSubjects matrix, e.g.
%  envCorrelation_z. 
%
%  Returns the parameter estimate from the contrast , beta
%  and the p-value associated with this under permutation

meanConn = get_mean_connectivity(netmat);


% Use NBS as a hacky shortcut to run the permuted GLM
GLM      = setup_glm(meanConn, Settings, iContrast);
testStat = NBSglm(GLM);

% parameter estimate and p-value
[beta, p] = estimate(testStat, GLM);

end



function meanConn = get_mean_connectivity(netmat)
% input dimensions
[nNodes, checkMe, nSubjects] = size(netmat);
assert(checkMe == nNodes,             ...
       [mfilename ':BadNetmatInput'], ...
	   'Please input 3D netmat. \n');
   
triUpperInd = logical(triu(ones(nNodes), 1));

% extract mean connectivity for each subject
for iSubject = nSubjects:-1:1,
	tmp                = netmat(:,:,iSubject);
	meanConn(iSubject) = mean(tmp(triUpperInd));
end%for

end%get_mean_connectivity

function GLM = setup_glm(m, Settings, iContrast)
%SETUP_GLM creates structure in appropriate size

GLM.y        = m(:);
GLM.perms    = 5000;
GLM.X        = Settings.GroupLevel.designMatrix; 
GLM.contrast = Settings.GroupLevel.contrasts(iContrast,:);
GLM.test     = 'ttest';
end%setup_glm

function [b, p] = estimate(testStat, GLM)

% compute the parameter estimate
b = GLM.contrast * (GLM.X \ GLM.y);

% compute p-value based on permutations of test statistic
p = mean(testStat(1) <= testStat(2:end));
end%estimate