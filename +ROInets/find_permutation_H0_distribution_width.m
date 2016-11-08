function sigma = find_permutation_H0_distribution_width(envData, nPerms, Regularize, transform)
%FIND_PERMUTATION_H0_DISTRIBUTION_WIDTH use random shuffles to build a null
%
% SIGMA = FIND_PERMUTATION_H0_DISTRIBUTION_WIDTH(ENVDATA, NPERMS, REGULARIZE, TRANSFORM)
%   finds the standard deviation of z-transformed correlations for null
%   data sharing the same temporal smoothness as input ENVDATA. NPERMS
%   complete datasets of the same size as ENVDATA are generated to build
%   this empirical null. 
%
%   If TRANSORM is true, a logarithmic transform is used to improve the
%   match between normally-generated variables and the input data. 
%
% SEE ALSO: FIND_EMPIRICAL_H0_DISTRIBUTION. 


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
%	$Revision: 368 $
%	$LastChangedDate: 2014-12-13 19:05:24 +0000 (Sat, 13 Dec 2014) $
%	Contact: giles.colclough 'at' magd.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 07-Nov-2013 12:43:45
if ~nPerms,
    % we will assume normality
    [sigma.z,         ...
     sigma.z_partial, ...
     sigma.z_partial_reg] = deal([]);
    return
end%if

if nargin < 3 || isempty(Regularize), 
    Regularize.do = false; 
end%if 
if Regularize.do,
    assert(isfield(Regularize, 'rho') && Regularize.rho > 0, ...
           [mfilename ':BadRegularizationParam'],            ...
           'Regularization parameter must be non-zero\n');
end%if
if nargin < 4
    transform = false;
end%if

nTrials     = size(envData, 3);
nPermsTotal = nTrials * nPerms;

for iPerm = nPermsTotal:-1:1,
    iTrial = mod(iPerm-1, nTrials)+1;
    surrogateData                                      = make_surrogate_data(envData(:,:,iTrial), ...
                                                                             transform);
    [corr(:,iPerm), pcorr(:,iPerm), pcorrReg(:,iPerm)] = compute_correlations(surrogateData, ...
                                                                              Regularize);
end%permutations

sigma.z             = std(ROInets.Fisher_r_to_z(corr(:)));
sigma.z_partial     = std(ROInets.Fisher_r_to_z(pcorr(:)));
sigma.z_partial_reg = std(ROInets.Fisher_r_to_z(pcorrReg(:)));

end%find_H0_permutation_distribution_width

function newData = make_surrogate_data(envData, transform)
%MAKE_SURROGATE_DATA generates one round of surrogates of ENVDATA
nModes = ROInets.rows(envData);

% make more normal
if transform, envData = log(envData); end

% phase randomisation
for iMode = nModes:-1:1,
    newData(:,iMode) = ROInets.generate_phase_surrogates(envData(iMode,:)',false,0); % Note transposition
end%for

% undo transform
if transform, newData = e.^newData; end
end%make_surrogate_data

function [r, pr, prR] = compute_correlations(data, Regularize)
%COMPUTE_CORRELATIONS generates correlation, partial correlation and
%  regularized partial correlation
%

upperTriInds = triu(true(ROInets.cols(data)),1);

rCov = cov(data);
cr   = corrcov(rCov, 1);
cpr  = ROInets.convert_precision_to_pcorr(ROInets.cholinv(rCov));
if Regularize.do,
    cprR = ROInets.convert_precision_to_pcorr(ROInets.dp_glasso(rCov, [], Regularize.rho));
else
    cprR = cpr;
end%if

r   = cr(upperTriInds);
pr  = cpr(upperTriInds);
prR = cprR(upperTriInds);
end%compute_correlations
