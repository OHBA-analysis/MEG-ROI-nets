function [nComponents, adjacency, pVals] = network_based_statistic_group_level(netmat, Settings, iContrast, threshold, alpha)
%network_based_statistic_group_level
%
% Runs a post-hoc NBS (Zalesky, 2012) on one netmat. 
%
% NETWORK_BASED_STATISTIC_GROUP_LEVEL(NETMAT, SETTINGS, CONTRAST_NO, THRESH, ALPHA) runs on
%   the nNodes x nNodes x nSubjects netmat (e.g. envCorrelation_z) using
%   Settings structure, choosing contrast CONTRAST_NO from the group-level
%   design. 
%
%   THRESH is the cut-off for the Z-stat on each edge. Try 2?
%   ALPHA is the FWER significance level. 
%   
%
% Sorry it's so clunky.  
%
% You need the NBS toolbox to run this: https://sites.google.com/site/bctnet/comparison/nbs


% Setup and run GLM
GLM      = setup_glm(netmat, Settings, iContrast);
testStat = NBSglm(GLM);

% show histogram of z-stats
plot_z_and_threshold(testStat(1,:), threshold);

% Setup and run stats
STATS = setup_stats(testStat, ROInets.rows(netmat), threshold, alpha);
[nComponents, sigComponents, pVals] = NBSstats(STATS);
for iC = nComponents:-1:1,
	tmp               = full(sigComponents{iC});
	adjacency(:,:,iC) = tmp + tmp.';
end%for

end%network_based_statistic_group_level









function GLM = setup_glm(netmat, Settings, iContrast)
%SETUP_GLM creates structure in appropriate size

GLM.y        = reformat_netmat(netmat);
GLM.perms    = 5000;
GLM.X        = Settings.GroupLevel.designMatrix;
GLM.contrast = Settings.GroupLevel.contrasts(iContrast,:);
GLM.test     = 'ttest';
end%setup_glm

function STATS = setup_stats(testStat, nModes, threshold, alpha)
%SETUP_STATS creates STATS parameters tructure for NBS
STATS.thresh    = threshold; % threshold for forming connected components
STATS.alpha     = alpha;     % significance threshold for FWER
STATS.N         = nModes;
STATS.test_stat = testStat;
STATS.size      = 'intensity';  % 'extent' or 'intensity' - measures either size of component, or total connectivity in component.
end%setup_stats

function y = reformat_netmat(data)
% turn netmats into a single row, taking only upper triangular part

[nNodes, checkMe, nSubjects] = size(data);
assert(checkMe == nNodes,             ...
       [mfilename ':BadNetmatInput'], ...
	   'Please input 3D netmat. \n');
   
triUpperInd = logical(triu(ones(nNodes), 1));

for iS = nSubjects:-1:1,
	tmp     = data(:,:,iS);
	y(iS,:) = tmp(triUpperInd);
end%for

end%reformat_netmat

function plot_z_and_threshold(z, threshold)
figure('Name', 'Z-stats input to NBS', 'Color', 'w');
hist(z, 30);
vline(threshold);
xlabel('Z');
end%plot