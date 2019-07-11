function [O, Q, R, W] = householder_orthogonalise(A)
%HOUSEHOLDER_ORTHOGONALISE	orthogonalisation using householder method
% O = HOUSEHOLDER_ORTHOGONALISE(A) Computes orthogonal set of vectors O 
%    from input A. Vectors must form columns of A. 
%
% [O, Q, R, W] = HOUSEHOLDER_ORTHOGONALISE(A) Computes orthonormal vectors
%    in Q, and R such that A = Q*R, and W such that O = A * W;
%
% See also: ROInets.symmetric_orthogonalise

%	Copyright 2013 OHBA
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
%	$Revision: 239 $
%	$LastChangedDate: 2014-08-15 14:58:49 +0100 (Fri, 15 Aug 2014) $
%	Contact: giles.colclough@gmail.com
%	Originally written on: GLNXA64 by Giles Colclough, 06-Nov-2013 13:35:17

[Q, R] = qr(A,0);

O = Q * (diag(diag(R)));

if nargout > 3,
	W = inv(R) * diag(diag(R)); %#ok<MINV> inverse acceptable for upper triangular matrix
end%if
