function correlationMats = do_group_level_glm(correlationMats, Settings)
%DO_GROUP_LEVEL_GLM run group comparison 
% 


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


%	$LastChangedBy: giles.colclough@gmail.com $
%	$Revision: 231 $
%	$LastChangedDate: 2014-08-07 20:53:06 +0100 (Thu, 07 Aug 2014) $
%	Contact: giles.colclough@gmail.com
%	Originally written on: MACI64 by Giles Colclough, 10-Apr-2014 13:33:21

for iFreq = Settings.nFreqBands:-1:1,
if strcmpi(Settings.paradigm, 'rest'),
    % univariate edge testing in turn for correlation, partial correlation
    % and regularized partials.
    [T, p, corrp, COPE] = ROInets.univariate_edge_test(correlationMats{iFreq}.envCorrelation_z, ...
                                                       Settings.GroupLevel.designMatrix,        ...
                                                       Settings.GroupLevel.contrasts);
    for iContrast = size(p,3):-1:1,
        h(:,:,iContrast) = ROInets.false_discovery_rate(ROInets.p_to_z_two_tailed(p(:,:,iContrast)), ...
                                                        Settings.FDRalpha);
    end%for
    
    groupLevel.correlation.T    = T;
    groupLevel.correlation.p    = p;
    groupLevel.correlation.FWEp = corrp;
    groupLevel.correlation.FDRh = h;
    groupLevel.correlation.COPE = COPE;
    
    [T, p, corrp, COPE] = ROInets.univariate_edge_test(correlationMats{iFreq}.envPartialCorrelation_z, ...
                                                 Settings.GroupLevel.designMatrix,        ...
                                                 Settings.GroupLevel.contrasts);
    for iContrast = size(p,3):-1:1,
        h(:,:,iContrast) = ROInets.false_discovery_rate(ROInets.p_to_z_two_tailed(p(:,:,iContrast)), ...
                                                        Settings.FDRalpha);
    end%for
    
    groupLevel.partialCorrelation.T    = T;
    groupLevel.partialCorrelation.p    = p;
    groupLevel.partialCorrelation.FWEp = corrp;
    groupLevel.partialCorrelation.FDRh = h;
    groupLevel.partialCorrelation.COPE = COPE;
    
    [T, p, corrp, COPE] = ROInets.univariate_edge_test(correlationMats{iFreq}.envPartialCorrelationRegularized_z, ...
                                                 Settings.GroupLevel.designMatrix,        ...
                                                 Settings.GroupLevel.contrasts);
    for iContrast = size(p,3):-1:1,
        h(:,:,iContrast) = ROInets.false_discovery_rate(ROInets.p_to_z_two_tailed(p(:,:,iContrast)), ...
                                                        Settings.FDRalpha);
    end%for
    
    groupLevel.partialCorrelationRegularized.T    = T;
    groupLevel.partialCorrelationRegularized.p    = p;
    groupLevel.partialCorrelationRegularized.FWEp = corrp;
    groupLevel.partialCorrelationRegularized.FDRh = h;
    groupLevel.partialCorrelationRegularized.COPE = COPE;
    
    % add in design for reference
    groupLevel.designMatrix = Settings.GroupLevel.designMatrix;
    groupLevel.contrasts    = Settings.GroupLevel.contrasts;
    
    
elseif strcmpi(Settings.paradigm, 'task'),
    % we need to run a separate GLM for each first level contrast. 
    % And for each group-level contrast.
    % What a mare.
    for iContrast = length(Settings.FirstLevel.contrasts):-1:1,
		% we use parameter estimates from the level below for each of
		% correlation, partial correlation and regularised partial
		% correlation
		if isfield(correlationMats{iFreq}, 'subjectLevel'),
			COPE = correlationMats{iFreq}.subjectLevel(iContrast).cope;
		elseif isfield(correlationMats{iFreq}, 'firstLevel'),
			COPE = correlationMats{iFreq}.firstLevel(iContrast).cope;
		else
			error([mfilename ':WhereIsTheData'], ...
				  'Expected input to have either first level or subject level results. \n');
		end%if
		
		
		%-- correlation
        [T, p, corrp, COPE] = ROInets.univariate_edge_test(COPE.correlation,                 ...
                                                           Settings.GroupLevel.designMatrix, ...
													       Settings.GroupLevel.contrasts);
        for iConG = size(p,3):-1:1,
            h(:,:,iConG) = ROInets.false_discovery_rate(ROInets.p_to_z_two_tailed(p(:,:,iConG)), ...
                                                        Settings.FDRalpha);
        end%for

        groupLevel(iContrast).correlation.T    = T;
        groupLevel(iContrast).correlation.p    = p;
        groupLevel(iContrast).correlation.FWEp = corrp;
        groupLevel(iContrast).correlation.FDRh = h;
        groupLevel(iContrast).correlation.COPE = COPE;

		%-- partial correlation
        [T, p, corrp, COPE] = ROInets.univariate_edge_test(COPE.partialCorrelation,          ...
                                                           Settings.GroupLevel.designMatrix, ...
                                                           Settings.GroupLevel.contrasts);
        for iConG = size(p,3):-1:1,
            h(:,:,iConG) = ROInets.false_discovery_rate(ROInets.p_to_z_two_tailed(p(:,:,iConG)), ...
                                                        Settings.FDRalpha);
        end%for

        groupLevel(iContrast).partialCorrelation.T    = T;
        groupLevel(iContrast).partialCorrelation.p    = p;
        groupLevel(iContrast).partialCorrelation.FWEp = corrp;
        groupLevel(iContrast).partialCorrelation.FDRh = h;
        groupLevel(iContrast).partialCorrelation.COPE = COPE;

		%-- regularised partial correlation
        if Settings.Regularize.do,
            [T, p, corrp, COPE] = ROInets.univariate_edge_test(COPE.partialCorrelationRegularized, ...
                                                               Settings.GroupLevel.designMatrix,   ...
                                                               Settings.GroupLevel.contrasts);
            for iConG = size(p,3):-1:1,
                h(:,:,iConG) = ROInets.false_discovery_rate(ROInets.p_to_z_two_tailed(p(:,:,iConG)), ...
                                                            Settings.FDRalpha);
            end%for

            groupLevel(iContrast).partialCorrelationRegularized.T    = T;
            groupLevel(iContrast).partialCorrelationRegularized.p    = p;
            groupLevel(iContrast).partialCorrelationRegularized.FWEp = corrp;
            groupLevel(iContrast).partialCorrelationRegularized.FDRh = h;
            groupLevel(iContrast).partialCorrelationRegularized.COPE = COPE;
        end%if
        % add in first levels for reference
        groupLevel(iContrast).firstLevelContrast        = Settings.SubjectLevel.contrasts{iContrast};
        groupLevel(iContrast).firstLevelConditionLabels = Settings.SubjectLevel.conditionLabel;
        
        % add in design for reference
        groupLevel(iContrast).groupDesignMatrix = Settings.GroupLevel.designMatrix;
        groupLevel(iContrast).groupContrasts    = Settings.GroupLevel.contrasts;
    end%for
    
else
    error([mfilename 'BadParadigm'], ...
          'Unrecognised paradigm %s. \n', Settings.paradigm);
end%if

correlationMats{iFreq}.groupLevel = groupLevel;
end%for loop over frequencies
