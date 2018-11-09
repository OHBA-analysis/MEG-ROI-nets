function [nodeData, voxelWeightings] = get_node_tcs(voxelData,spatialBasis,timeCourseGenMethod,ROIname,ROIlabels) 
% GET_NODE_TCS extracts ROI time-courses
%
% INPUTS
% - voxelData: (nVoxels x nSamples) matrix *OR* a filename for a .mat file
%   containing one variable with a (nVoxels x nSamples) matrix *OR* a SPM MEEG
%   object
% - spatialBasis: (nVoxels x nParcels) logical matrix of parcel membership OR
%   (nVoxels x nParcels) matrix with parcel membership weights
% - timeCourseGenMethod: one of {'PCA','peakVoxel','spatialBasis','mean'}.
% - ROI_name: if voxelData is an MEEG object, write new montage with this
%   string within the name
% - ROIlabels: if voxelData is an MEEG object, label channels with this cell
%   array (must contain nParcels elements)
%
% OUTPUTS
% - nodeData: (nParcels x nSamples) matrix of node timecourses if voxelData
%   was a matrix, *or* an SPM MEEG object with a new online montage containing
%   nParcels channels
% - voxelWeightings: relative weighting of voxels over the brain, used to
%   construct the time-course for each ROI.
%
% USAGE PATTERNS
%
% There are broadly two usage choices
%
% - voxelData could be
%   - a matrix, in which case the output (nodeData) is a matrix
%   - an MEEG object, in which case the output (nodeData) is an MEEG object
%     with a new online montage
% - spatialBasis could be
%   - parcel membership, in which case spatialBasis should be a logical matrix
%     consisting of 1s and 0s and timecourseGenMethod can be
%     {'PCA','peakVoxel','mean'}
%   - A spatial basis set (e.g. from group ICA), in which case spatialBasis
%     has voxel weights for each parcel that can be weighted and overlapping,
%     and timecourseGenMethod should be set to 'spatialBasis'
%
% TIMECOURSE METHODS
% 
% timeCourseGenMethod chooses how the ROI time-course is created. It can be:
%   - 'PCA'   - take 1st PC of voxels
%   - 'peakVoxel' - voxel with the maximum variance
%   - 'mean' - NOT RECOMMENDED. Mean timecourse of voxels. This method is not
%     sign-invariant which makes it invalid for most beamformer data.
%   - 'spatialBasis' - Infer timecourses of a spatial basis set that is
%     (nVoxels x nSamples). The ROI time-course for each spatial map is the
%     1st PC from all voxels, weighted by the spatial map. If the parcellation
%     is unweighted and non-overlapping, 'spatialBasis' will give the same
%     result as 'PCA' except with a dfiferent normalization
%
% SPATIAL BASIS INPUT SPECIFICATION
% 
% If 'spatialBasis' is a parcel membership matrix, then it should be a logical
% array (nVoxels x nParcels) which identifies the voxels in each parcel. Each
% column identifies the voxels making up each parcel. True entries indicate
% membership.
% 
% If 'spatialBasis' is a spatial basis set (e.g. from group ICA) then each
% spatial map (held in columns) is a whole-brain map - each map can be non-
% binary and maps may overlap. Note that the spatial basis should be
% orthogonal, or nearly so. This is because parcel data are computed one
% parcel at a time, without multiple regression, which is valid for Data = SB
% * NodeData + e only if SB are orthogonal.
%
% SPATIAL BASIS ALGORITHM
%
% First scale all maps so they have a positive peak of height 1. Then for each
% map in the spatial basis:
% 
%     1. Variance-normalise the clean voxel time-series. (This is now turned
%        off.) Weight by the clean spatial map.
% 
%     2. Perform a PCA. Extract the coefficients of the first PC to represent
%        the node time-course
% 
%     3. Re-weight the node time-course to match the variance in the ROI.
%    
%
% Copyright 2014-2017 OHBA, FMRIB This program is free software: you can
% redirstribute it and/or modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation, either version 3 of
% the License, or (at your option) any later version.
%	
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
% more details.
%	
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.

%% find representative time courses for each parcel

if ischar(voxelData) && strcmp(voxelData(end-3:end), '.mat') && exist(voxelData, 'file'), 
    tmp       = load(voxelData);
    ff        = fieldnames(tmp);
    voxelData = tmp.(ff{1});
    clear tmp;
    goodSamples = true(1,ROInets.cols(voxelData));

elseif ischar(voxelData),
    error([mfilename ':FileInputError'],               ...
          ['Did not recognise the filename input. \n', ...
           '  Ensure the file exists and has a .mat extension. ']);

elseif isa(voxelData,'meeg')
    goodSamples = good_samples(voxelData);
else
    % data passed in
    goodSamples = true(1,ROInets.cols(voxelData));
end%if

nParcels = ROInets.cols(spatialBasis);

% sort out names for ROIs
if nargin < 5 || isempty(ROIlabels),
  ROIlabels = arrayfun(@(x) sprintf('ROI %d',x),1:nParcels,'UniformOutput',false);
else
	assert(iscell(ROIlabels) && length(ROIlabels) == nParcels, ...
		[mfilename ':BadROIlabels'],                         ...
		'ROIlabels should be a cell array matching the number of ROIs in the parcellation. \n');
end%if

if nargin < 4 || isempty(ROIname) 
  ROIname = 'Parcellated';
end


ft_progress('init', 'text', '');

if isa(voxelData,'meeg')
  assert(voxelData.nchannels==size(spatialBasis,1),sprintf('Parcellation has %d voxels, but data has %d channels. Is the correct montage selected?',size(spatialBasis,1),voxelData.nchannels));
else
  assert(size(voxelData,1)==size(spatialBasis,1),sprintf('Parcellation has %d voxels, but data has %d rows. Is the data correctly oriented?',size(spatialBasis,1),size(voxelData,1)));
end


switch lower(timeCourseGenMethod)
    case 'pca'
        if any(spatialBasis(:)~=0 & spatialBasis(:)~=1),
            warning([mfilename ':NonBinaryParcelMask'],    ...
                    ['Input parcellation is not binary. ', ...
                     'Parcel weights will be removed. \n']);
        end%if
        spatialBasis = logical(spatialBasis);    
        % check that each voxel is only a member of one parcel
        assert(~any(ROInets.row_sum(spatialBasis) > 1), ...
               [mfilename ':MultipleParcelOccupancy'], ...
               'Each voxel can be a member of at most one parcel. \n');

        % demean and variance normalise each voxel - 20 May 2014 Don't
        % variance normalise as MEG data are very smooth, and power is a
        % good indication of true sources
        if isa(voxelData,'meeg')
            temporalSTD = max(sqrt(osl_source_variance(voxelData)), eps);
        else
            temporalSTD = max(std(voxelData, [], 2), eps);
        end
        
        % pre-allocate PCA weightings for each parcel
        voxelWeightings = zeros(size(spatialBasis));
        
        % perform PCA on each parcel and select 1st PC scores to represent
        % parcel
        for iParcel = nParcels:-1:1,
            progress = nParcels - iParcel + 1;
            ft_progress(progress / nParcels, ...
                        [mfilename ...
                         ':    Finding PCA time course for ROI %d out of %d'], ...
                        iParcel, nParcels);
                
            thisMask = spatialBasis(:, iParcel);
            if any(thisMask), % non-zero
                parcelData = voxelData(find(thisMask),:,:); %#ok Can't use logical indexing
                parcelData = parcelData(:,goodSamples);
                parcelData = ROInets.demean(parcelData, 2);
                which_nan = isnan(parcelData);
                parcelData(which_nan) = 0;
                [U, S, V]  = ROInets.fast_svds(parcelData, 1);
                PCAscores  = S * V';
                
                % restore sign and scaling of parcel time-series
                % U indicates the weight with which each voxel in the
                % parcel contributes to the 1st PC
                TSsign          = sign(mean(U));
                relVoxelWeights = abs(U) ./ sum(abs(U)); % normalise the linear combination
                % weight the temporal STDs from the ROI by the proportion used in 1st PC
                TSscale         = dot(relVoxelWeights, temporalSTD(thisMask)); 
                nodeTS          = TSsign .*                               ...
                                  (TSscale / max(std(PCAscores), eps)) .* ... 
                                  PCAscores;
                      
                % return the linear operator which is applied to the data
                % to retrieve the nodeTS
                voxelWeightings(thisMask, iParcel) = TSsign .* ...
                                                     (TSscale / max(std(PCAscores), eps)) ...
                                                     .* U';
                
            else
                warning([mfilename ':EmptySpatialComponentMask'],          ...
                        ['%s: When calculating ROI time-courses, ',        ...
                         'an empty spatial component mask was found for ', ...
                         'component %d. \n',                               ...
                         'The ROI will have a flat zero time-course. \n',  ...
                         'Check this does not cause further problems ',    ...
                         'with the analysis. \n'],                         ...
                         mfilename, iParcel);
                     
                %nodeTS = zeros(1, ROInets.cols(voxelData));
		nodeTS = zeros(1, size(voxelData,2)*size(voxelData,3));
            end%if
            
            nodeData(iParcel,:) = nodeTS;
        end%for
        
        clear parcelData 

    case 'peakvoxel'
        if any(spatialBasis(:)~=0 & spatialBasis(:)~=1),
            warning([mfilename ':NonBinaryParcelMask'],    ...
                    ['Input parcellation is not binary. ', ...
                     'It will be binarised. \n']);
        end%if
        spatialBasis = logical(spatialBasis);
        
        % find rms power in each voxel
        if isa(voxelData,'meeg')
            voxelPower = osl_source_variance(voxelData); % I'm going to use variance instead, because it's fast to compute
        else
            voxelPower = sqrt(ROInets.row_sum(voxelData.^2) ./ ...
                              ROInets.cols(voxelData));
        end
        
        % pre-allocate weightings for each parcel
        voxelWeightings = zeros(size(spatialBasis));
                      
        % take peak voxel in each parcel
        for iParcel = nParcels:-1:1,
            progress = nParcels - iParcel + 1;
            ft_progress(progress / nParcels, ...
                        [mfilename ...
                         ':    Finding peak voxel time course for ROI %d out of %d'], ...
                        iParcel, nParcels);
            
            thisMask = spatialBasis(:, iParcel);
            
            if any(thisMask), % non-zero
                % find index of voxel with max power
                thisParcPower            = voxelPower;
                thisParcPower(~thisMask) = 0;
                [~, maxPowerInd]         = max(thisParcPower);
                
                % select voxel timecourse
                nodeData(iParcel,:) = voxelData(maxPowerInd,:);
                
                % save which voxel was used
                voxelWeightings(maxPowerInd, iParcel) = 1;
                
            else
                warning([mfilename ':EmptySpatialComponentMask'],          ...
                        ['%s: When calculating ROI time-courses, ',        ...
                         'an empty spatial component mask was found for ', ...
                         'component %d. \n',                               ...
                         'The ROI will have a flat zero time-course. \n',  ...
                         'Check this does not cause further problems ',    ...
                         'with the analysis. \n'],                         ...
                         mfilename, iParcel);
                     
                nodeData(iParcel,:) = zeros(1, ROInets.cols(voxelData));
            end%if
        end%loop over parcels
        
        clear parcelData
        
    case 'spatialbasis'
        % scale group maps so all have a positive peak of height 1
        % in case there is a very noisy outlier, choose the sign from the
        % top 5% of magnitudes
        top5pcInd = abs(spatialBasis) >=                        ...
                         repmat(prctile(abs(spatialBasis), 95), ...
                                [ROInets.rows(spatialBasis), 1]);
        for iParcel = nParcels:-1:1,
            mapSign(iParcel) = sign(mean(...
                              spatialBasis(top5pcInd(:,iParcel), iParcel)));
        end%for
        scaledSpatialMaps = ROInets.scale_cols(spatialBasis, ...
                                   mapSign ./                ...
                                   max(max(abs(spatialBasis), [], 1), eps));
        
        % estimate temporal-STD for normalisation
        if isa(voxelData,'meeg')
            temporalSTD = max(sqrt(osl_source_variance(voxelData)), eps);
        else
            temporalSTD = max(std(voxelData, [], 2), eps);
        end

        voxelWeightings = zeros(size(spatialBasis)); % Preallocate the nvoxels x nparcels weight matrix      
        
        % find time-course for each spatial basis map
        for iParcel = nParcels:-1:1, % allocate memory on the fly
            progress = nParcels - iParcel + 1;
            ft_progress(progress / nParcels, ...
                        [' ' mfilename ...
                         ':    Finding spatial basis time course for ROI %d out of %d'], ...
                        iParcel, nParcels);
            
            % extract the spatial map of interest
            thisMap     = scaledSpatialMaps(:, iParcel);
            parcelMask  = logical(thisMap);
                        
            % variance-normalise all voxels to remove influence of
            % outliers. - remove this step 20 May 2014 for MEG as data are
            % smooth and little risk of high-power outliers. Also, power is
            % a good indicator of sensible signal. 
            % Weight all voxels by the spatial map in question
            % AB - apply the mask first then weight, to reduce memory use
            weightedTS  = voxelData(find(parcelMask),:,:); %#ok Can't use logical indexing
            weightedTS  = reshape(weightedTS,[size(weightedTS,1),size(weightedTS,2)*size(weightedTS,3)]); % reshape to handle trials
            weightedTS  = weightedTS(:,goodSamples);
            weightedTS  = ROInets.scale_rows(weightedTS, thisMap(parcelMask));

            % perform svd and take scores of 1st PC as the node time-series
            % U is nVoxels by nComponents - the basis transformation
            % S*V holds nComponents by time sets of PCA scores - the 
            % timeseries data in the new basis
            [U, S, V]   = ROInets.fast_svds(weightedTS, 1);
            clear weightedTS
            
            PCAscores   = S * V';
            maskThresh  = 0.5; % 0.5 is a decent arbitrary threshold chosen by Steve Smith and MJ after playing with various maps.
            thisMask    = thisMap(parcelMask) > maskThresh;   
            
            if any(thisMask), % the mask is non-zero
                % U is the basis by which voxels in the mask are weighted
                % to form the scores of the 1st PC
                relativeWeighting = abs(U(thisMask)) ./ ...
                                    sum(abs(U(thisMask)));
                
                TSsign  = sign(mean(U(thisMask)));
                TSscale = dot(relativeWeighting, temporalSTD(thisMask));       
                nodeTS  = TSsign .*                               ...
                          (TSscale / max(std(PCAscores), eps)) .* ...      
                          PCAscores;
                      
                % for Mark: this is the linear operator which is applied to
                % the voxel data to get nodeTS.
                voxelWeightings(parcelMask,iParcel) = TSsign .* ...
                                             (TSscale / max(std(PCAscores), eps)) ...
                                             .* (U' .* thisMap(parcelMask)');
                
            else
                warning([mfilename ':EmptySpatialComponentMask'],          ...
                        ['%s: When calculating ROI time-courses, ',        ...
                         'an empty spatial component mask was found for ', ...
                         'component %d. \n',                               ...
                         'The ROI will have a flat zero time-course. \n',  ...
                         'Check this does not cause further problems ',    ...
                         'with the analysis. \n'],                         ...
                         mfilename, iParcel);
                     
                nodeTS = zeros(1, ROInets.cols(weightedTS));
                voxelWeightings(~thisMask, iParcel) = zeros(length(thisMask), 1);
            end%if
            
            nodeData(iParcel, :) = nodeTS;
            
        end%loop over parcels
        
    case 'mean'
      % WARNING - Sign flipping will affect the mean timecourse
      % Use with caution - generally PCA is preferable to mean
      
      % Parcel timecourse will be the weighted average of voxels belonging to that parcel
      voxelWeightings = bsxfun(@rdivide,spatialBasis,sum(spatialBasis));

      if ~isa(voxelData, 'meeg')
        nodeData = voxelWeightings'*voxelData;
      end

    otherwise
        error([mfilename ':UnrecognisedTimeCourseMethod'],            ...
              ['Unrecognised method for finding ROI time-course. \n', ...
               'Expected ''PCA'', ''spatialBasis'', or ''mean''. \n']);
end%switch

ft_progress('close');

if isa(voxelData, 'meeg')
  % convert output into an online montage
  % Assume if ROInets is being used with SPM12 this is happening through OSL
  % i.e. add_montage() is present
  currentMontage      = voxelData.montage('getmontage');
  if isempty(currentMontage)
    name = ROIname;
  else
    name = [ROIname ' - ' currentMontage.name];
  end

  nodeData = add_montage(voxelData,voxelWeightings.',name,ROIlabels);

end
