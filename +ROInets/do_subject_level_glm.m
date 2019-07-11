function correlationMats = do_subject_level_glm(correlationMats, Settings)
%DO_SUBJECT_LEVEL_GLM run subject level analysis
%
% M = ROInets.DO_SUBJECT_LEVEL_GLM(M, Settings) uses subject design
%   matrix Settings.SubjectLevel.subjectDesign to perform a subject-level analysis, and append these results
%   as a subject level to the set of correlation matrices M. 
%
%   The subject design matrix should be an nSessions by nSubjects matrix
%   with entries 1 or zeros to indicate allocation of sessions to subjects.
%   
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

if ~strcmpi(Settings.paradigm, 'task'),
	% no need to do anything: only 1 session assumed for RS analysis
	return
end%if
    
% quick dimension check
nSessions = size(correlationMats{1}.firstLevel(1).cope.correlation,3);

if ~isfield(Settings, 'SubjectLevel') || isempty(Settings.SubjectLevel) || isempty(Settings.SubjectLevel.subjectDesign),
	subjectDesign = eye(nSessions);
else
	subjectDesign = Settings.SubjectLevel.subjectDesign;
end%if

assert(ROInets.rows(subjectDesign) == nSessions, ...
	   [mfilename ':dimensionMismatch'], ...
	   'subjectDesign should have as many rows as there are sessions to analyse. \n');
    
for iFreq = Settings.nFreqBands:-1:1,    
    % we need to run a separate GLM for each first level contrast. 
    for iContrast = length(Settings.FirstLevel.contrasts):-1:1, 
		
		% we use parameter estimates from the level below for each of
		% correlation, partial correlation and regularised partial
		% correlation
		if isfield(correlationMats{iFreq}, 'firstLevel'),
			COPE = correlationMats{iFreq}.firstLevel(iContrast).cope;
		else
			error([mfilename ':WhereIsTheData'], ...
				  'Expected input to have first level results. \n');
		end%if
		
		% find any subjects which have nan entries. They will need to be
		% removed from the data and design matrix. 
		% It's a safe assumption that NaNs encode subjects without
		% sufficient information in trials to estimate parameters 
		% They will therefore be nan over all edges, and for each of
		% correlation, partialcorr and the regularised version
		nanSessions  = any(isnan(COPE.correlation),1);
		goodSessions = logical(squeeze(~nanSessions(1,1,:)));
		cleanDesign  = subjectDesign(goodSessions,:);
		
		
		%-- correlation
		beta = cleanDesign \ ROInets.get_edges(COPE.correlation(:,:,goodSessions))';
		subjectLevel(iContrast).cope.correlation = ROInets.unvectorize(beta');

		%-- partial correlation
        beta = cleanDesign \ ROInets.get_edges(COPE.partialCorrelation(:,:,goodSessions))';
		subjectLevel(iContrast).cope.partialCorrelation = ROInets.unvectorize(beta');

		%-- regularised partial correlation
        if Settings.Regularize.do,
            beta = cleanDesign \ ROInets.get_edges(COPE.partialCorrelationRegularized(:,:,goodSessions))';
			subjectLevel(iContrast).cope.partialCorrelationRegularized = ROInets.unvectorize(beta');
        end%if
        % add in first levels for reference
        subjectLevel(iContrast).firstLevelContrast        = Settings.SubjectLevel.contrasts{iContrast};
        subjectLevel(iContrast).firstLevelConditionLabels = Settings.SubjectLevel.conditionLabel;
        
        % add in design for reference
        subjectLevel(iContrast).designMatrix = subjectDesign;
    end%for
    


correlationMats{iFreq}.subjectLevel = subjectLevel;
end%for loop over frequencies
