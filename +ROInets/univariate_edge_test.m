function [T, p, corrp, COPE] = univariate_edge_test(netmats, designMatrix, contrasts, standardise)
%UNIVARIATE_EDGE_TEST permutation test for group-level stats on networks
%
% [T, P, CORRP] = UNIVARIATE_EDGE_TEST(NETMATS, DESIGN, CONTRAST) performs
%    univariate testing for significance on each edge of a network matrix. 
%    Testing is performed using FSL's randomise, and uses 5000 permutations
%    of the group labels to perform nonparametric inference. 
%
%    Pass in NETMATS, which are symmetric network matrices with subjects in
%    the third dimension. The diagonals will be ignored. The DESIGN matrix
%    should have as many rows as subjects. The CONTRAST matrix should have
%    as many columns as the design matrix has columns. 
%
% [T, P, CORRP] = UNIVARIATE_EDGE_TEST(..., STANDARDISE) demeans and
%    variance normalises the design matrix, if TRUE. 
%
% [T, P, CORRP, COPE] = ... also returns the contrast of parameter
%   estimates from the regression. 
%
%   T provides single-edge T-stats; P the uncorrected p-Values; and CORRP the
%   FWE-corrected p-values. To use the weaker FDR correction, see
%   ROINETS.FALSE_DISCOVERY_RATE. 
%
% This code demeans over subjects at the top level, so it is not good for
% testing the whole-group mean effect.

%	Copyright 2015 OHBA, FMRIB
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


%	$LastChangedBy: GilesColclough $
%	$Revision: 763 $
%	$LastChangedDate: 2015-10-21 11:52:19 +0100 (Wed, 21 Oct 2015) $
%	Contact: giles.colclough@gmail.com
%	Originally written on: MACI64 by Giles Colclough, 25-Sep-2014 15:20:18


%% Input checking
[nNodes, checkMe, nSessions] = size(netmats);
assert(checkMe == nNodes, ...
       [mfilename ':NonsquareInput'], ...
       'Input netmats must be square, with subjects in the third dimension. \n');
   
[checkMe, nEVs] = size(designMatrix);
assert(checkMe == nSessions, ...
       [mfilename ':BadDesign'], ...
       'Design matrix must have as many rows as subjects. \n');
   
assert(ROInets.cols(contrasts) == nEVs, ...
       [mfilename ':BadContrasts'], ...
       'Contrasts must have as many columns as EVs in the design matrix. \n');

if nargin < 4 || ~exist('standardise', 'var'),
    standardise = false;
else
    assert(islogical(standardise), ...
           [mfilename ':BadStandardise'], ...
           'Standardise input must be a logical value. \n');
end%if
   
edges = ROInets.get_edges(netmats); % note this assumes symmetry

[Ttmp, ptmp, corrptmp, COPEtmp] = ROInets.perform_glm_with_randomise(edges,        ...
														             designMatrix, ...
                                                                     contrasts,    ...
                                                                     standardise);

% convert back to symmetric matrices
T     = ROInets.unvectorize(Ttmp);
p     = ROInets.unvectorize(1 - ptmp);
corrp = ROInets.unvectorize(1 - corrptmp);
COPE  = ROInets.unvectorize(COPEtmp);

end%univariate_edge_tests
% [EOF]
