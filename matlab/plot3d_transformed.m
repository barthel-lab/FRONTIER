% MRI_Coords_Plot_offline.m - plots the coordinates for MRI biopsy data
% Script assumes metadata excel sheet has already been loaded to workspace
% as a table named VUmc

% Add include directory to path
addpath('include')
addpath('include/spatialmath')

% File pointers
nodfp = '../results/meta/FRONTIER.nodes.transformed.csv'
segfp = '../results/meta/FRONTIER.seg.transformed.csv'
niimapfp = '../results/meta/FRONTIER.nii.csv'
outmovbase = '../results/mov/'
vidfp = '../results/mov/test'
hvoxcortfp = '/usr/local/fsl/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz'
hvoxsubfp = '/usr/local/fsl/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr0-1mm.nii.gz'

% Load node positions and segments connecting nodes
noddat = readtable(nodfp)
segdat = readtable(segfp)

% Identify unique patients in dataset
pts = unique(noddat.Patient)

% Loop ALL patients
for pti = 1:length(pts)
    
    % Print figure for each patient
    figure
    hold on
    
    plotCortex()
    plotSubcortex()
    plotLatVent()

    % Set patient to plot volumes for
    patient = char(pts(pti))
    
    % Initialize node and segment data for selected patient
    nodes = noddat(find(strcmp(noddat.Patient,patient)),:)
    seg = segdat(find(strcmp(segdat.Patient,patient)),:)
    
    % Plot tumor masks for a given patient
    plotTumor(patient)
    
    % Export movie file
    exportMov([patient,'.mask'],patient)
    
    % Plot phylo 3D for given patient
    plotPhylo3D(nodes,seg)
     
    % Export movie file
    exportMov([patient,'.mask.phylo'],patient)
end

% Manually set a number of patients to plot together
selected_pts = {'VUmc-01', 'VUmc-05', 'Vumc-17'} %, 'VUmc-07', 'VUmc-11', 'Vumc-15'}
fn = ''

figure
hold on
    
plotCortex()
plotSubcortex()
plotLatVent()

% Loop across manually selected patients and plot together
for pti = 1:length(selected_pts)

    % Set patient to plot volumes for
    patient = char(selected_pts(pti))
    
    fn = [fn, '.', patient]
    
    % Initialize node and segment data for selected patient
    nodes = noddat(find(strcmp(noddat.Patient,patient)),:)
    seg = segdat(find(strcmp(segdat.Patient,patient)),:)
    
    % Plot tumor masks for a given patient
    plotTumor(patient)
    
    % Plot phylo 3D for given patient
    plotPhylo3D(nodes,seg)
    
    % Export movie file
    exportMov(['agg.phylo',fn], '')

end


% Plot all patients together
figure
hold on
    
plotCortex()
plotSubcortex()
plotLatVent()

% Loop ALL patients
for pti = 1:length(pts)

    % Set patient to plot volumes for
    patient = char(pts(pti))
    
    % Initialize node and segment data for selected patient
    nodes = noddat(find(strcmp(noddat.Patient,patient)),:)
    seg = segdat(find(strcmp(segdat.Patient,patient)),:)
    
    % Plot phylo 3D for given patient
    plotPhylo3D(nodes,seg)
end

% Export movie file
exportMov('agg.ALL','')

% Draw cortex + white matter
figure
plotCortex()
plotWM()
exportMov('brain_wm','Cortex + WM')

% Draw cortex + subcortical nuclei
figure
plotCortex()
plotSubcortex()
exportMov('cortex_subcortex','Cortex + Subcortex')
plotLatVent()
exportMov('cortex_subcortex_latvent','Cortex + Subcortex + Ventricles')
plotSubcortex2()
exportMov('cortex_subcortex_latvent_full','Cortex + Full Subcortex + Ventricles')

% END %