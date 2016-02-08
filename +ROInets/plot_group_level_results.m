function plot_group_level_results(correlationMats, Settings, netType, correctionType, ROInames, ROIorder)
%PLOT_GROUP_LEVEL_RESULTS
%
% PLOT_GROUP_LEVEL_RESULTS(CORRELATIONMATS, SETTINGS, NETTYPE, CORRECTIONTYPE, 
%                          ROINAMES, ROIORDER)
%   plots results from group-level network analysis. 
%
%   Takes in CORRELATIONMATS, output from osl_network_analysis. Plots
%   heatmaps of significance using multiple comparisons CORRECTIONTYPE
%   for network analysis method NETTYPE. Plots will be labelled using cell
%   array ROINAMES, with optional permutation of ROIs ROIORDER. 
%
%   NETTYPE: - 'correlation'
%			 - 'partialCorrelation'
%			 - 'partialCorrelationRegularized'
%
%   CORRECTIONTYPE - 'T'    - plot T-stats
%				   - 'p'    - plot uncorrected p-values
%				   - 'FDRh' - threshold uncorrected p-values at FDR = 0.05
%				   - 'FWEp' - plot family-wise error corrected p-values
%
%   Be aware: p-values are displayed as -log10(p). 
%        e.g. 0.05  -> 1.3
%			  0.01  -> 2
%		      0.001 -> 3
%
%
% Example usage:
% 
% fmri_d100reduced_labels();
% ROInets.plot_grou_level_results(correlationmats, Settings, 'correlation', 'FWEp', LABELS(NEW_ORDER), NEW_ORDER);

%% parse inputs
if ~iscell(correlationMats),
	error([mfilename ':NonCellInput'], ...
		  'Expecting a cell array of structures as first input. \n');
end%if
isTask = isfield(correlationMats{1}, 'firstLevel') ...
	             && strcmpi(Settings.paradigm, 'task');
if ~isTask,
	error([mfilename ':NotTaskInput'], ...
		  'Expecting correlationMats to be output of task analysis. \n');
end%if	

nFreqs                    = Settings.nFreqBands;
nFreqsCheck               = length(correlationMats);
nFirstLevelContrasts      = length(Settnings.SubjectLevel.contrasts);
nFirstLevelContrastsCheck = length(correlationMats{1}.firstLevel);
nGroupLevelContrasts      = ROInets.rows(Settings.GroupLevel.contrasts);
nGroupLevelContrastsCheck = size(correlationMats{1}.groupLevel(1).correlation.T,3);

assert(nFreqsCheck               == nFreqs               && ...
	   nFirstLevelContrastsCheck == nFirstLevelContrasts && ...
	   nGroupLevelContrastsCheck == nGroupLevelContrasts,   ...
	   [mfilename ':MatSettingsMismatch'],                  ...
	   'The Settings file and CorrelationMats file do not seem to match up. \n');

nROIs  = ROInets.rows(correlationMats{1}.firstLevel(1).cope.correlation);
validatestring(netType, {'correlation', 'partialCorrelation', ...
					     'partialCorrelationRegularized'},    ...
	           mfilename, 'netType', 3);
validatestring(correctionType, {'T', 'p', 'FDRh', 'FWEp'}, ...
	           mfilename, 'correctionType', 4);
if nargin >= 5 && ~isempty(ROInames),
	assert((iscell(ROInames) && length(ROInames) == nROIs));
else
	ROInames = 1:nROIs;
end%if
if nargin >= 6 && ~isempty(ROIorder),
	assert(length(ROIorder) == nROIs);
else
	ROIorder = 1:nROIs;
end%if

%% Let's plot away
FONTSIZE = 15;
for iFreq = 1:nFreqs,
	for iConFirst = 1:nFirstLevelContrasts,
			matBase = correlationMats{iFreq}.groupLevel(iConFirst).(netType);
		for iConGroup = 1:nGroupLevelContrasts,
			CAXIS = [0 3];
			switch correctionType
				case 'T'
			        plotMat = matBase.(correctionType)(ROIorder,ROIorder,iConGroup);
					CAXIS = [0 15];
					
				case {'p', 'FWEp'}
					plotMat = -log10(matBase.(correctionType)(ROIorder,ROIorder,iConGroup));
					
				case 'FDRh'
					% threshold original p-values by FDR-corrected bound
					h = matBase.(correctionType)(ROIorder,ROIorder,iConGroup);
					p = matBase.p(ROIorder,ROIorder,iConGroup);
					threshP       = p;
					threshP(p~=h) = 1;
					plotMat       = -log10(threshP);
				otherwise
					error([mfilename ':BadCorrectionType'], ...
						  'Unrecognised correction type. \n');
			end
			plotTitle = sprintf('Freq band %d, First level contrast %d, Group level contrast %d, %s %s results', ...
				                iFreq, iConFirst, iConGroup, netType, correctionType);
		    figure('Name', plotTitle, 'Color', 'w');
			imagesc(plotMat(ROIorder,ROIorder));
			colormap(bluewhitered);
			set(gca, 'XTick', 1:nROIs, 'XTickLabel', ROInames);
			set(h,...
				'FontName', 'Helvetica', ...
				'FontSize', FONTSIZE, ...
				'Box', 'on', ...
				'YGrid', 'off', ...
				'XGrid', 'off', ...
				'TickDir', 'in', ...
				'TickLength', [0.005 0.005], ...
				'XMinorTick', 'off', ...
				'YMinorTick', 'off', ...
				'XColor', [0.3 0.3 0.3], ...
				'YColor', [0.3 0.3 0.3], ...
				'LineWidth', 2);			
			colorbar;
		end%for
	end%for
end%for


end%plot_froup_level_results
% [EOF]