function [L, d, rho, W] = closest_orthogonal_matrix(A)
%CLOSEST_ORTHOGONAL_MATRIX  Computes closest orthogonal matrix
%
% L = CLOSEST_ORTHOGONAL_MATRIX(A) returns orthogonal matrix L which
%   is closest to A, as measured by the Frobenius norm of (L-A), such that
%   transpose(L) * L is non-negative diagonal. 
%
% [L, D, RHO] = CLOSEST_ORTHOGONAL_MATRIX(A) also returns the scaling
%   factors for each column of A, D, and the final converged square
%   distance between L and A, RHO. 

%	References:
%	R. Evans,
%	"http://empslocal.ex.ac.uk/people/staff/reverson/uploads/Site/procrustes.pdf"
%   "Orthogonal, but not Orthonormal, Procrustes Problems", 1997
%	
%   See also: ROInets.symmetric_orthogonalise, svd. 

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

%	$LastChangedBy$
%	$Revision$
%	$LastChangedDate$
%	Contact: giles.colclough@gmail.com
%	Originally written on: MACI64 by Giles Colclough, 07-Aug-2014 13:27:52
%	Optimised by JH on Jul 19th 2016

% settings
MAX_ITER  = 2e2;
TOLERANCE = max(1, max(size(A)) * ROInets.fast_svds(A,1)) * eps(class(A)); % the tolerance used for matrices within matlab
DEBUG     = false;

% tests for convergence
reldiff     = @(a,b) 2*abs(a-b) / (abs(a)+abs(b));
absdiff     = @(a,b) abs(a-b); % previously used - JH
convergence = @(rho, prevRho) reldiff(rho,prevRho) <= TOLERANCE;

% declare memory and initialise
iter = 0;
A_b  = conj(A);
d    = sqrt(dot( A, A_b, 1 )); % alternative: ones( 1, ROInets.cols(A) )
rho  = NaN(MAX_ITER, 1);

while iter < MAX_ITER,
    iter = iter + 1;
    
    % find orthonormal polar factor
    try
        V = ROInets.symmetric_orthogonalise(ROInets.scale_cols( A, d ));
    catch ME
        % use suppression of some variables such that rank is reduced as a
        % stopping criterion
        if iter > 1 && strcmp(ME.identifier, 'ROInets:RankError'),
            % stop the process and use that last L computed.
            rank_reduction_warning();
            break
        else
            throw(ME); % Re-throw the original error
        end%if
    end%try
    
    % minimise rho = |A - V D|^2 w.r.t d
    d = dot( A_b, V, 1 );
    
    % new best orthogonal estimate
    L = ROInets.scale_cols( V, d );
    
    % compute error term
    E         = A - L;
    rho(iter) = sqrt(sum(dot( E, conj(E), 1 ))); % use sqrt(trace) instead - JH
    
    if DEBUG,
        fprintf('%s: iter %4d \t rho %g \n', mfilename, iter, rho(iter));  %#ok<UNRCH>
    end
    
    if iter > 1 && convergence(rho(iter), rho(iter - 1)),
        break
    end%if
end%convergence loop

if nargout > 3,
    % output linear operator - run once more
    [~, ~, ~, W] = ROInets.symmetric_orthogonalise(ROInets.scale_cols( A, d ));
    W = diag(d) * W * diag(d);
end%if

% tidy vector
rho(isnan(rho)) = [];

if isequal(iter, MAX_ITER),
    warning([mfilename ':MaxIterationsHit'],                              ...
        ['%s: hit max iterations: %d. ',                              ...
        'You may not have found the optimal orthogonal matrix. \n'], ...
        mfilename, MAX_ITER);
end%if

if DEBUG,
    figure('Name', 'Convergence path', 'Color', 'w');                      %#ok<UNRCH>
    semilogy(2:(length(rho)-1), rho(2:end-1) - rho(end), 'b+');
    xlabel('iteration');
    ylabel('\rho_i - \rho_{\inf}');
end%if

end%closest_orthogonal_matrix

function [] = rank_reduction_warning()
% This warning is produced when the orthogonalisation process runs out of
% rank. 
% This can happen when two vectors are very close to each other but
% different in magnitude. The alignement process can suppress the smaller
% vector down to zero, as there is little residual left. This, in some
% circumstances, can reduce the rank below the dimensionality of the data.
%
% As L at all iterations of the algorithm are orthogonal matrices, and as
% after one iteration the magnitudes do become close to the original set,
% it seems like this is a fair stopping criterion.
warning([mfilename ':RankReducedTooFar'],                                 ...
		['%s: Orthogonalisation optimisation stopped when rank reduced ', ...
         'to limit. \n'],                                                 ...
		mfilename);
end%rank_reduction_warning
% [EOF]














% % % 
% % % % test some methods for computation within the iteration loop
% % % 
% % % %% Setup
% % % nNodes   = 40;
% % % nSamples = 200*600;
% % % nRuns    = 20;
% % % 
% % % % make covariance matrix
% % % sigmaOff = 0.2 * rand(nNodes);
% % % sigmaOff = sigmaOff' * sigmaOff;
% % % sigma = sigmaOff - diag(diag(sigmaOff)) + eye(nNodes);
% % % 
% % % var = 1 + 0.1 * randn(nNodes, 1);
% % % sigma = diag(sqrt(var)) * sigma * diag(sqrt(var));
% % % 
% % % % make mean matrix
% % % mu = 4 + randn(1, nNodes);
% % % 
% % % % simulate data stream
% % % A = mvnrnd(mu, sigma, nSamples);
% % % A = ROInets.demean(A);
% % % 
% % % V = ROInets.symmetric_orthogonalise(A);
% % % %% A finding a d vector
% % % fprintf('finding d\n');
% % % 
% % % % 1. diag
% % % diag_time = zeros(1,nRuns);
% % % d = zeros(nNodes,1);
% % % 
% % % for iRun = 1:nRuns,
% % %     tic
% % %     d = diag(A' * V);
% % %     diag_time(iRun) = toc;
% % % end
% % % 
% % % fprintf('Diag time: %g\n', sum(diag_time));
% % % fprintf('max d: %g\n\n', max(d));
% % % 
% % % % 2. vector norms
% % % loop_time = zeros(1,nRuns);
% % % for iRun = 1:nRuns,
% % %     tic
% % %     for k = nNodes:-1:1,
% % %         d(k) = A(:,k)' * V(:,k);
% % %     end
% % %     loop_time(iRun) = toc;
% % % end
% % % fprintf('Loop time: %g\n', sum(loop_time));
% % % fprintf('max d: %g\n\n', max(d));
% % % 
% % % %% B finding rho
% % % fprintf('finding rho\n');
% % % d = diag(A' * V);
% % % rho = NaN;
% % % % % 1. all in one trace
% % % % AIO_time = zeros(1,nRuns);
% % % % for iRun = 1:nRuns,
% % % %     tic
% % % %     rho = sum(diag( (A - V * diag(d))' * (A - V * diag(d)) )); % inline trace
% % % %     AIO_time(iRun) = toc;
% % % % end
% % % % fprintf('All in one trace time: %g\n', sum(AIO_time));
% % % 
% % % % 2. pre-allocate trace
% % % AIO_time = zeros(1,nRuns);
% % % for iRun = 1:nRuns,
% % %     tic
% % %     M = A - V * diag(d);
% % %     rho = sum(diag( M' * M )); % inline trace
% % %     AIO_time(iRun) = toc;
% % % end
% % % fprintf('Pre-allocate trace time: %g\n', sum(AIO_time));
% % % fprintf('rho: %g\n\n', rho);
% % % 
% % % % on the spot
% % % OTS_time = zeros(1,nRuns);
% % % for iRun = 1:nRuns,
% % %     tic
% % %     M = A - V * diag(d);
% % %     rho = sum(sum(M.^2)); % inline trace
% % %     OTS_time(iRun) = toc;
% % % end
% % % fprintf('On-the-spot time: %g\n', sum(OTS_time));
% % % fprintf('rho: %g\n\n', rho);
% % % 
% % % % 3. norm
% % % norm_time = zeros(1,nRuns);
% % % for iRun = 1:nRuns,
% % %     tic
% % %     rho = norm( (A - V * diag(d)), 'fro' )^2; 
% % %     norm_time(iRun) = toc;
% % % end
% % % fprintf('Norm time: %g\n', sum(norm_time));
% % % fprintf('rho: %g\n\n', rho);
% % % 
% % % % % 4. svd
% % % % svd_time = zeros(1,nRuns);
% % % % for iRun = 1:nRuns,
% % % %     tic
% % % %     rho = sum(svd(A - V * diag(d)).^2); 
% % % %     svd_time(iRun) = toc;
% % % % end
% % % % fprintf('SVD time: %g\n', sum(svd_time));
% % % 
% % % % 5. loop
% % % loop_time = zeros(1,nRuns);
% % % for iRun = 1:nRuns,
% % %     tic
% % %     for k = nNodes:-1:1,
% % %         rho_(k) = sum( (A(:,k) - V(:,k) * d(k)).^2 );
% % %     end
% % %     rho = sum(rho_); 
% % %     loop_time(iRun) = toc;
% % % end
% % % fprintf('Loop time: %g\n', sum(loop_time));
% % % fprintf('rho: %g\n\n', rho);
% % % 
% % % %% Check to see if this way is better
% % % [L, d, rhoFinal] = ROInets.closest_orthogonal_matrix(A);
% % % display(rho - rhoFinal(end));
