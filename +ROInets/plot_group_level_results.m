function plot_group_level_results(correlationMats, Settings, netType, correctionType, ROInames, ROIorder, firstLevelConsToPlot, groupLevelConsToPlot)
%PLOT_GROUP_LEVEL_RESULTS
%
% PLOT_GROUP_LEVEL_RESULTS(CORRELATIONMATS, SETTINGS, NETTYPE, CORRECTIONTYPE, 
%                          ROINAMES, ROIORDER, FIRSTLEVELCONSTOPLOT, 
%						   GROUPLEVELCONSTOPLOT)
%   plots results from group-level network analysis. 
%
%   Takes in CORRELATIONMATS, output from run_network_analysis. Plots
%   heatmaps of significance using multiple comparisons CORRECTIONTYPE
%   for network analysis method NETTYPE. Plots will be labelled using cell
%   array ROINAMES, with optional permutation of ROIs ROIORDER. Chooses
%   which FIRSTLEVELCONSTOPLOT and GROUPLEVELCONSTOPLOT. 
%
%   The last four options can be omitted or set to [] to receive defaults
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

nFreqs                    = Settings.nFreqBands;
nFreqsCheck               = length(correlationMats);
nGroupLevelContrasts      = ROInets.rows(Settings.GroupLevel.contrasts);
nGroupLevelContrastsCheck = size(correlationMats{1}.groupLevel(1).correlation.T,3);

assert(nFreqsCheck               == nFreqs               && ...
	   nGroupLevelContrastsCheck == nGroupLevelContrasts,   ...
	   [mfilename ':MatSettingsMismatch'],                  ...
	   'The Settings file and CorrelationMats file do not seem to match up. \n');

if isTask,
   nROIs  = ROInets.rows(correlationMats{1}.firstLevel(1).cope.correlation);
else
   nROIs = ROInets.rows(correlationMats{1}.correlation);
end%if

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
if nargin >= 8 && ~isempty(groupLevelConsToPlot),
	assert(all(ismember(groupLevelConsToPlot, 1:nGroupLevelContrasts)));
else
	groupLevelConsToPlot = 1:nGroupLevelContrasts;
end%if


	FONTSIZE = 12; % MP edit from 15
	

if isTask,
	nFirstLevelContrasts      = length(Settings.SubjectLevel.contrasts);
	nFirstLevelContrastsCheck = length(correlationMats{1}.firstLevel);
	assert(nFirstLevelContrastsCheck == nFirstLevelContrasts, ...
	       [mfilename ':FirstLevelMismatch'],                 ...
	       'The Settings file and CorrelationMats file do not seem to match up. \n');
	if nargin >= 7 && ~isempty(firstLevelConsToPlot),
		assert(all(ismember(firstLevelConsToPlot, 1:nFirstLevelContrasts)));
	else
		firstLevelConsToPlot = 1:nFirstLevelContrasts;
	end%if

	%% Let's plot away
	for iFreq = 1:nFreqs,
		for iConFirst = firstLevelConsToPlot,
				matBase = correlationMats{iFreq}.groupLevel(iConFirst).(netType);
			for iConGroup = groupLevelConsToPlot,
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
				set(gca, 'YTick', 1:nROIs, 'YTickLabel', ROInames);
				set(gca,...
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
				caxis(CAXIS);
				colormap(bluewhitered);
				axis square
				colorbar;
			end%for
		end%for
	end%for

else % not task
	%% Let's plot away
	for iFreq = 1:nFreqs,
		matBase = correlationMats{iFreq}.groupLevel.(netType);
		for iConGroup = groupLevelConsToPlot,
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
			plotTitle = sprintf('Freq band %d, Group level contrast %d, %s %s results', ...
								iFreq, iConGroup, netType, correctionType);
			figure('Name', plotTitle, 'Color', 'w');
			imagesc(plotMat(ROIorder,ROIorder));
			set(gca, 'YTick', 1:nROIs, 'YTickLabel', ROInames);
            set(gca, 'XTick', 1:nROIs, 'XTickLabel', ROInames);
			set(gca,...
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
			caxis(CAXIS);
			colormap(bluewhitered);
			axis square
			colorbar;
		end%for
	end%for
end%if istask

end%plot_froup_level_results
% [EOF]
