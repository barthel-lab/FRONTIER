% MRI_Coords_Plot_offline.m - plots the coordinates for MRI biopsy data
% Script assumes metadata excel sheet has already been loaded to workspace
% as a table named VUmc

% Add matlab to path
addpath('matlab')

% Load node positions and segments connecting nodes
noddat = readtable('FRONTIER.nodes.csv')
segdat = readtable('FRONTIER.seg.csv')

% Identify unique patients in dataset
pts = unique(noddat.Patient)

for pti = 1:length(pts)

    % Set patient to plot volumes for
    patient = char(pts(pti))
    
    % Initialize node and segment data for selected patient
    snoddat = noddat(find(strcmp(noddat.Patient,patient)),:)
    ssegdat = segdat(find(strcmp(segdat.Patient,patient)),:)
    
    % Nii file directory
    filedir = 'nii';
    
    % Load nii filemap
    niifiles = readtable('FRONTIER.nii.csv')
    sniifiles = niifiles(find(strcmp(niifiles.Patient,patient)),:)
    
    % Prepare nii file paths
    t1gfile = ''
    if(strlength(sniifiles.t1g) > 0)
        t1gfile = char(strcat(filedir, '/', sniifiles.t1g));
    end
    flairfile = ''
    if(strlength(sniifiles.flair) > 0)
        flairfile = char(strcat(filedir, '/', sniifiles.flair))
    end
    maskfile = char(strcat(filedir, '/', sniifiles.mask))
    samplefile = char(strcat(filedir, '/', sniifiles.samples))
    
    % Extract sample coordinates
    % [~,s] = extractSamples(samplefile);
    
    % Get number of samples
    j = height(snoddat(string(snoddat.Tip)=="TRUE",:))
    
    % Faster alternative to extract sample coordinates
    s = table2array(snoddat(1:j,{'X','Y','Z'}))
    
    % Extract brain mask
    brainmask = extractCoords(maskfile);
    
    % Calculate tumor distances, corresponding points on tumor surface,
    % centroid, radius, corresponding point on tumor surface, and volume
    [BrainDistance, ~, BrainCentroid, BrainRadius, ~, BrainVol] = calcDist(s, brainmask);
    
    % Calculate distances between points and centroid
    BrainCentroidDistance = calcEucDist(BrainCentroid, s)
    
    % Expand radius and volume for output
    BrainRadius = repmat(BrainRadius,j,1)
    BrainVol = repmat(BrainVol,j,1)
    
    % Extract T1G volume
    if(strlength(t1gfile) > 0)
        t1gmask = extractCoords(t1gfile);
        [T1GDistance, ~, T1GCentroid, T1GRadius, ~, T1GVol] = calcDist(s, t1gmask);
        T1GCentroidDistance = calcEucDist(T1GCentroid, s)
        T1GRadius = repmat(T1GRadius,j,1)
        T1GVol = repmat(T1GVol,j,1)
    else
        T1GDistance = repmat(missing,j,1);
        T1GCentroidDistance = repmat(missing,j,1);
        T1GRadius = repmat(missing,j,1);
        T1GVol = repmat(missing,j,1);
    end
    
    % Extract FLAIR volume
    if(strlength(flairfile) > 0)
        flairmask = extractCoords(flairfile);
        [FLAIRDistance, ~, FLAIRCentroid, FLAIRRadius, ~, FLAIRVol] = calcDist(s, flairmask);
        FLAIRCentroidDistance = calcEucDist(FLAIRCentroid, s)
        FLAIRRadius = repmat(FLAIRRadius,j,1)
        FLAIRVol = repmat(FLAIRVol,j,1)
    else
        FLAIRDistance = repmat(missing,j,1);
        FLAIRCentroidDistance = repmat(missing,j,1);
        FLAIRRadius = repmat(missing,j,1);
        FLAIRVol = repmat(missing,j,1);
    end
    
    % Extract biopsy order
    Biopsy = snoddat.Biopsy(1:j)
    
    % Repeat patient
    Patient = repmat(patient,j,1)
    
    % Create output table
    out = table(Patient, Biopsy, BrainDistance, BrainCentroidDistance, BrainRadius, BrainVol, T1GDistance, T1GCentroidDistance, T1GRadius, T1GVol, FLAIRDistance, FLAIRCentroidDistance, FLAIRRadius, FLAIRVol)
    
    if(exist('res'))
        res = vertcat(res,out)
    else
        res = out
    end
    
end

writetable(res, 'FRONTIER.distances.csv')