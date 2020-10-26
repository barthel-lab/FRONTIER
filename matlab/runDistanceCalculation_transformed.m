% MRI_Coords_Plot_offline.m - plots the coordinates for MRI biopsy data
% Script assumes metadata excel sheet has already been loaded to workspace
% as a table named VUmc

addpath('include')
addpath('include/spatialmath')

% File pointers
nodfp = '../results/meta/FRONTIER.nodes.transformed.csv'
segfp = '../results/meta/FRONTIER.seg.transformed.csv'

% Load node positions and segments connecting nodes
nodes = readtable(nodfp)
seg = readtable(segfp)

% Identify unique patients in dataset
pts = unique(nodes.Patient)

% Extract volumes
cort = plotCortex()
wmat = plotWM()
[thal, caud, puta] = plotSubcortex()
vent = plotLatVent()
[pall, hipp, amyg, accu] = plotSubcortex2()

% Calculate distances between points and subcortal structures
s = table2array(nodes(string(nodes.Tip)=="TRUE",{'X','Y','Z'}))
[dist_cort, ~, cent_cort, ~, ~, ~] = calcDist(s, cort);
[dist_wmat, ~, ~, ~, ~, ~] = calcDist(s, wmat);
[dist_thal, ~, ~, ~, ~, ~] = calcDist(s, thal);
[dist_caud, ~, ~, ~, ~, ~] = calcDist(s, caud);
[dist_puta, ~, ~, ~, ~, ~] = calcDist(s, puta);
[dist_vent, ~, ~, ~, ~, ~] = calcDist(s, vent);
[dist_pall, ~, ~, ~, ~, ~] = calcDist(s, pall);
[dist_hipp, ~, ~, ~, ~, ~] = calcDist(s, hipp);
[dist_amyg, ~, ~, ~, ~, ~] = calcDist(s, amyg);
[dist_accu, ~, ~, ~, ~, ~] = calcDist(s, accu);

% Calculate distances between points and cortical centroid
dist_cent_cort = calcEucDist(cent_cort, s)
    
%nodes(1:j,{'Patient', 'Biopsy'}), 
res = horzcat(nodes(string(nodes.Tip)=="TRUE",{'Patient', 'Biopsy', 'X', 'Y', 'Z'}), table(dist_cort, dist_wmat, dist_thal, dist_caud, dist_puta, dist_vent, dist_pall, dist_hipp, dist_amyg, dist_accu))

writetable(res, '../results/meta/FRONTIER.distances.transformed.csv')

% Backup nodes
%nodes_all = nodes

% Verify points
figure
hold on
plotCortex()
plotSubcortex()
%plotWM()

% Select points inside WM
sres = res%(res.dist_wmat>0,:)

% Plot nodes in 3d space
for i = 1:size(sres,1)
    x = sres.X(i)
    y = sres.Y(i)
    z = sres.Z(i)
    plot3(x, y, z, 'o', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black', 'MarkerSize', 1, 'LineWidth', 2.5)
end

% END %