function nodeData = remove_source_leakage(nodeDataOrig, protocol, output_meeg)
%REMOVE_SOURCE_LEAKAGE correct ROI time-courses for source leakage
%
% NODEDATA = REMOVE_SOURCE_LEAKAGE(NODEDATAORIG, PROTOCOL)
%   produces orthogonalised node time-courses NODEDATA from 
%   uncorrected node time-courses NODEDATAORIG
%
%   PROTOCOL is a string to switch between various all-to-all 
%   orthogonalisation methods for source-spread correction. It can be:
%     'none'          - No correction applied. 
%     'symmetric'     - Apply orthogonalisation on the parcel time-courses.
%                       This produces orthonormal parcel time-courses
%                       which are as close as possible to the original
%                       time-courses.
%     'closest'       - Apply orthogonalisation on the parcel time-courses.
%                       Start as for the symmetric method, then converge to
%                       a (not orthonormal) orthogonal matrix which is as
%                       close as possible to the original time-courses. 
%     'householder'   - Orthogonalise using a more numerically stable
%                       alternative to the Gram-Schmidt process, dealing
%                       with ROI time-courses in a random order. 
%
    

%	Copyright 2014 OHBA
%	This program is free software: you can redirstribute it and/or modify
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
%	$Revision: 235 $
%	$LastChangedDate: 2014-08-07 22:09:15 +0100 (Thu, 07 Aug 2014) $
%	Contact: giles.colclough 'at' magd.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough and Stephen Smith.



% catch cases where the rank of the node data is lower than the number of
% nodes being orthogonalised
rankErrorMessage = ['The ROI time-course matrix is not full rank. \n',    ...
                    '    This prevents you from using an all-to-all ',    ...
                    'orthogonalisation method. \n',                       ...
                    '    Your data have rank %d, and you are looking ',   ...
                    'at %d ROIs. \n',                                     ...
                    '    You could try reducing the number of ROIs, or ', ...
                    'using an alternative orthogonalisation method. \n'];

if nargin < 3 || isempty(output_meeg) 
  output_meeg = true;
end

switch protocol
    case 'none'
        % no orthogonalisation applied to parcel time-courses
        nodeData=nodeDataOrig;
        
    case 'closest'
        % finds closest orthogonal set of vectors by applying symmetric 
        % orthogonalisation then iterating to find closest orthogonal matrix
        orthogFunctionHandle = @ ROInets.closest_orthogonal_matrix;
        
        nodeData = find_closest_orthogonal_matrix(nodeDataOrig,         ...
                                                  rankErrorMessage,     ...
                                                  orthogFunctionHandle, ...
												  protocol);
        
    case 'symmetric'
        % finds closest orthonormal matrix
        orthogFunctionHandle = @ ROInets.symmetric_orthogonalise;
        
        nodeData = find_closest_orthogonal_matrix(nodeDataOrig,         ...
                                                  rankErrorMessage,     ...
                                                  orthogFunctionHandle, ...
												  protocol);
                                               
        
    case 'householder'
        nodeData = find_orthogonal_matrix_by_householder_reflections(...
                                           nodeDataOrig, rankErrorMessage);
        
    otherwise
        error([mfilename ':UnrecognisedOrthMethod'], ...
              'Unrecognised parcel orthogonalisation protocol. \n');
end%switch
end%remove_source_leakage
% -------------------------------------------------------------------------
function nodeData = find_closest_orthogonal_matrix(nodeDataOrig,     ...
                                                   rankErrorMessage, ...
                                                   orthogFunction,   ...
												   protocol)
%FIND_CLOSEST_ORTHOGONAL_MATRIX wrapper on orthogonalisation functions
nParcels = ROInets.rows(nodeDataOrig);

try
    fprintf('    Orthogonalising...\n');
	if isa(nodeDataOrig, 'meeg'),
		[~, ~, ~, W] = orthogFunction(transpose(nodeDataOrig(:,:)));   
		nodeData     = add_montage(nodeDataOrig, W, protocol);
		
	else
		nodeData = transpose(orthogFunction(transpose(nodeDataOrig)));
	end%if
	
% catch rank deficiency in node data
catch ME
    if regexp(ME.identifier, 'RankError'),
        error([mfilename ':RankError'], ...
              rankErrorMessage,         ...
              rank(nodeDataOrig), nParcels);
    else
        rethrow(ME);
    end%if
end%try
end%find_closest_orthogonal_matrix
% -------------------------------------------------------------------------
function nodeData = find_orthogonal_matrix_by_householder_reflections(...
                                            nodeDataOrig, rankErrorMessage)
%FIND_ORTHOGONAL_MATRIX_BY_HOUSEHOLDER_REFLECTIONS wrapper on householder
% orthogonalisation

nParcels = ROInets.rows(nodeDataOrig);

isRankDeficient = rank(nodeDataOrig) < nParcels;
if isRankDeficient,
    error([mfilename ':RankError'], ...
          rankErrorMessage,         ...
          rank(nodeDataOrig), nParcels);
end%if

fprintf('    Orthogonalising...\n');

permutation                     = randperm(nParcels);
permutationInverse(permutation) = 1:nParcels;

if isa(nodeDataOrig, 'meeg') && output_meeg
	% bugger the permutation bit for this faff
	[~, ~, ~, W] = ROInets.householder_orthogonalise(nodeDataOrig(:,:).'); 
	nodeData     = add_montage(nodeDataOrig, W, 'householder');
else
	nodeData = ROInets.householder_orthogonalise(nodeDataOrig(permutation,:).').'; 
	nodeData = nodeData(permutationInverse,:);
end

end%find_orthogonal_matrix_by_householder_reflections


function D = add_montage(D, W, protocol)
%ADD_MONTAGE adds the mapping W' to object D.
% set up new montage by premultiplication
nMontages       = D.montage('getnumber');
currentMontage  = D.montage('getmontage');
newMontage      = currentMontage;
unit            = unique(D.units());
newMontage.tra  = W.' * currentMontage.tra;
newMontage.name = [protocol ' orthogonalised ' currentMontage.name];

% convert to spm object
D = D.montage('add', newMontage);
D = D.montage('switch', nMontages + 1);
D = D.units(:,unit{1});
D.save();
end%add_montage
% [EOF]