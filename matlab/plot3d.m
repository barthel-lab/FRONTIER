% MRI_Coords_Plot_offline.m - plots the coordinates for MRI biopsy data
% Script assumes metadata excel sheet has already been loaded to workspace
% as a table named VUmc

% Add include directory to path
addpath('include')

% File pointers
nodfp = '../results/meta/FRONTIER.nodes.transformed.csv'
segfp = '../results/meta/FRONTIER.seg.transformed.csv'
niimapfp = '../results/meta/FRONTIER.nii.csv'
outmovbase = '../results/mov/'

% Load node positions and segments connecting nodes
noddat = readtable(nodfp)
segdat = readtable(segfp)

% Identify unique patients in dataset
pts = unique(noddat.Patient)

for pti = 1:length(pts)

    % Set patient to plot volumes for
    patient = char(pts(pti))
    
    % Initialize node and segment data for selected patient
    snoddat = noddat(find(strcmp(noddat.Patient,patient)),:)
    ssegdat = segdat(find(strcmp(segdat.Patient,patient)),:)
    
    % Video file outputs
    vidf1 = ['mov/',patient,'-3dlineage'];
    vidf2 = ['mov/',patient,'-3dlineage.tumor'];
    vidf3 = ['mov/',patient,'-3dlineage.tumor.brain'];
    
    % Initialize figure
    % figure
    % hold on
    
    % Plot segments connecting nodes
    for i = 1:size(ssegdat,1)
        p1 = snoddat(snoddat.Node==ssegdat.Node1(i),:)
        p2 = snoddat(snoddat.Node==ssegdat.Node2(i),:)
        plot3([p1.X; p2.X], [p1.Y; p2.Y], [p1.Z; p2.Z], '--black', 'MarkerFaceColor', 'black', 'MarkerSize', 15, 'LineWidth', 1)
    end
    
    % Plot nodes in 3d space
    for i = 1:size(snoddat,1)
        x = snoddat.X(i)
        y = snoddat.Y(i)
        z = snoddat.Z(i)
        label = snoddat.Biopsy(i)
        plot3(x, y, z, 'o', 'MarkerFaceColor', findColor(snoddat.Subtype(i)), 'MarkerEdgeColor', 'black', 'MarkerSize', findSize(string(snoddat.Tip(i))), 'LineWidth', 2.5)
        if(strcmp(label,'NA') == 0)
            text(x+0.3, y+0.3, z+0.3, label)
        end
    end
    
    % Axis equalize
    axis equal
    
    % Change Axis Font for clarity; add labels
    set(gca,'linewidth',2)
    set(gca,'fontsize',14)
    set(get(gca,'XLabel'),'string','X (mm)')
    set(get(gca,'YLabel'),'string','Y (mm)')
    set(get(gca,'ZLabel'),'string','Z (mm)')
    
    % Set title
    set(get(gca,'title'),'string', patient)
    
    % Set background ('w' = white, 'none' = transparent)
    set(gcf,'color','w')
    
    % Save video (3d lineage)
    % OptionZ.FrameRate=5;OptionZ.Duration=15;OptionZ.Periodic=true; % Rendering settings
    % CaptureFigVid([-20,10;-110,35;-190,80;-290,10;-380,10],vidf1,OptionZ)
    
    % Nii file directory
    filedir = 'nii';
    
    % Load nii filemap
    niifiles = readtable(niimapfp)
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
    
    % Draw T1G volume
    if(strlength(t1gfile) > 0)
        tumor1 = extractCoords(t1gfile);
        tumor1.FaceColor = 'red'
        tumor1.EdgeColor = 'none';
        tumor1.FaceAlpha = 0.20;
    end
    
    % Draw FLAIR volume
    if(strlength(flairfile) > 0)
        tumor2 = extractCoords(flairfile);
        tumor2.FaceColor = [1 0.65 0];
        tumor2.EdgeColor = 'none';
        tumor2.FaceAlpha = 0.15;
    end
    
    % Save video (3d lineage + tumor volume)
    % CaptureFigVid([-20,10;-110,35;-190,80;-290,10;-380,10],vidf2,OptionZ)
    
    % Add brain overlay
    hold on
    brain_overlay = extractCoords(maskfile)
    brain_overlay.FaceColor = 'black';
    brain_overlay.EdgeColor = 'none';
    brain_overlay.FaceAlpha = 0.08;
    
    % Save video (3d lineage + tumor volume + brain volume)
    % CaptureFigVid([-20,10;-110,35;-190,80;-290,10;-380,10],vidf3,OptionZ)

end
