function D = add_montage(D, W, name,ROIlabels)
	%ADD_MONTAGE adds the mapping W' to object D.
	% set up new montage by premultiplication
	%
	% This function is partly developed and is currently unused - finishing implementation 
	% and using it in remove_source_leakage() and get_node_tcs() is on the todo list

	if nargin < 4 || isempty(ROIlabels) 
		ROIlabels = [];
	end
	
	nMontages       = D.montage('getnumber');
	currentMontage  = D.montage('getmontage');
	newMontage      = currentMontage;
	unit            = unique(D.units());
	keyboard
	newMontage.tra  = W.' * currentMontage.tra;
	newMontage.name = name;

	if ~isempty(ROIlabels)
		newMontage.labelnew = ROIlabels;
	end

	% convert to spm object
	D = D.montage('add', newMontage);
	D = D.montage('switch', nMontages + 1);
	D = D.units(:,unit{1});
