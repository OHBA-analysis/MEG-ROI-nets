function P = tikhonov_inverse(C, rho)
%TIKHONOV_INVERSE inverse covariance using Tikhonov regularisation
%
% Regularises the estimate of the inverse covariance by adding rho*I 
%
% P = TIKHONOV_INVERSE(C, RHO) returns the inverse covariance P estimated
%    from sample covariance matrix C, using regularisation parameter RHO. 


%	Copyright 2016 OHBA
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
%	$Revision: 214 $
%	$LastChangedDate: 2014-07-24 12:40:42 +0100 (Thu, 24 Jul 2014) $
%	Contact: giles.colclough@gmail.com
%	Originally written on: MACI64 by Giles Colclough, 14-Jul-2016

if nargin < 2 || ~exist('rho', 'var') || isempty(rho),
	rho = 0.1; % Steve's default
end

nModes = ROInets.rows(C);
P      = GC_cholinv(C + rho * eye(nModes));
