function c = pli(phase)
	% Take in phase timecourses and compute phase lag index
	%
	% INPUT
	% - phase - phase timecourse, n_signals x n_times 
	% 
	% OUTPUT
	% - c - PLI connectivity matrix
	%
	% EXAMPLE USAGE
	%
	% 	pli(phase)
	%
	% Romesh Abeysuriya 2017

	% Remove any NaNs (or Infs)
	clean = all(isfinite(phase),1);
	phase = phase(:,clean);

	assert(all(isfinite(phase(:)))); 
	plv_out = zeros(size(phase,1));
	[a,b] = meshgrid(1:size(phase,1));
	a = triu(a,1);
	b = triu(b,1);
	a = a(:);
	b = b(:);
	idx = find(a~=0 & b~=0);

	c = zeros(size(phase,1));
	phase_diffs = zeros(size(phase,2),length(idx));

	for l = 1:length(idx)
		phase_diffs(:,l) = phase(a(idx(l)),:)-phase(b(idx(l)),:);
	end

	tmp = abs(mean(sign(sin(phase_diffs)),1));

	for l = 1:length(idx)
		c(a(idx(l)),b(idx(l))) = tmp(l);
	end

	c = tril(c)+tril(c)';

		
