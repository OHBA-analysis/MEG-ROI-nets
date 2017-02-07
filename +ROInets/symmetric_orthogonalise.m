function [L, ignore, ignore2, W] = symmetric_orthogonalise(A, maintainMagnitudes)
%SYMMETRIC_ORTHOGONALISE closest orthogonal matrix
% 
% L = SYMMETRIC_ORTHOGONALISE(A) returns orthonormal matrix L which
%   is closest to A, as measured by the Frobenius norm of (L-A). 
%
%   The orthogonal matrix is constructed from a singular value decomposition
%   of A. 
%
% L = SYMMETRIC_ORTHOGONALISE(A, KEEP_MAGNITUDES) returns the orthogonal
%   matrix L, whose columns have the same magnitude as the respective
%   columns of A, and which is closest to A, as measured by the Frobenius
%   norm of (L-A), if KEEP_MAGNITUDES is TRUE. 
%
%   The orthogonal matrix is constructed from a singular value decomposition
%   of A. 
%
% [L, ~, ~, W] = SYMMETRIC_ORTHOGONALISE(...) also returns a weighting matrix
%   such that L = A * W;
%
%   See also: ROINETS.HOUSEHOLDER_ORTHOGONALISE, ORTH, SVD. 

% References: Naidu, A. R. "Centrality of Lowdin Orthogonalizations",
%   arXiv 1105.3571v1, May 2011. 
%   Available at: http://arxiv-web3.library.cornell.edu/pdf/1105.3571v1.pdf

%	Copyright 2013 OHBA
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
%	$Revision: 239 $
%	$LastChangedDate: 2014-08-15 14:58:49 +0100 (Fri, 15 Aug 2014) $
%	Contact: giles.colclough 'at' magd.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 31-Oct-2013 13:30:05

if nargin < 2 || ~exist('maintainMagnitudes', 'var'),
    maintainMagnitudes = false;
end

[ignore, ignore2] = deal([]); % to match up with other functions

if maintainMagnitudes,
    D = diag(sqrt(diag(A' * A)));
    
    
	
	if nargout > 1,
		% call function again
		[Lnorm, ~, ~, W] = ROInets.symmetric_orthogonalise(A * D, false);
		
		% scale result
		L = Lnorm * D;
		W = D * W * D;
		
	else
		% call function again
		Lnorm = ROInets.symmetric_orthogonalise(A * D, false);
		
		% scale result
		L = Lnorm * D;
	end%if
	
else
    [U, S, V] = svd(A, 'econ');

    if ~isempty(S),
        % we need to check that we have sufficient rank
        S   = diag(S);
        tol = max(size(A)) * S(1) * eps(class(A));
        r   = sum(S > tol);
        isFullRank = (r >= ROInets.cols(A));

        if isFullRank,
            % polar factors of A
            L = U * conj(transpose(V));
			
			if nargout > 1,
				W = V * diag(1.0 ./ S) * conj(transpose(V));
			end%if
        else % not enough degrees of freedom
        	error_message = ['The ROI time-course matrix is not full rank. \n',    ...
        	                    '    This prevents you from using an all-to-all ',    ...
        	                    'orthogonalisation method. \n',                       ...
        	                    '    Your data have rank %d, and you are looking ',   ...
        	                    'at %d ROIs. \n',                                     ...
        	                    '    You could try reducing the number of ROIs, or ', ...
        	                    'using an alternative orthogonalisation method. \n'];
        	error('ROInets:RankError',error_message,r,ROInets.cols(A));
        end%if        

    else
        L = [];
		W = [];
    end%if
end%if
end%symmetric_orthogonalise
% [EOF]