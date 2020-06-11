% plotFigCalculations.m - Create a figure showing how tumor sample
% distances are calculated for MRI data
% Input is a number indicating which patient is to be plotted

% After rotating to prefered view, print using
% print(gcf,'Example_Calculation-P1','-dpdf','-r600');

function plotFigCalculations(patient)

% Identify the tumor core and sample files for a given patient
if (class(patient) == 'double') % default class for numbers is double
    patient = num2str(patient);
end

fileDirectory = '/Users/anderk/Documents/FRONTIER/VUmc Patient MRI Coords/Nii_Files';

switch patient
    case "1"
        coreFilename = [fileDirectory,'/','P01T1G.nii.gz'];
        sampleFilename = [fileDirectory,'/','P01samples.nii.gz'];
    case "2"
        coreFilename = 'P02T1G.nii.gz';
        sampleFilename = 'P02samples.nii.gz';
    case "3"
        coreFilename = 'P03_new_volume.nii.gz';
        sampleFilename = 'P03samples.nii.gz';
    case "4"
        coreFilename = 'P04T1G.nii.gz';
        sampleFilename = 'P04samples.nii.gz';
    case "5"
        coreFilename = 'P05T1G.nii.gz';
        sampleFilename = 'P05samples.nii.gz';
    case "6"
        coreFilename = 'P06T1G.nii.gz';
        sampleFilename = 'P06samples.nii.gz';
    case "7"
        coreFilename = 'P07T1G.nii.gz';
        sampleFilename = 'P07samples.nii.gz';
    case "8"
        coreFilename = 'P08T1G.nii.gz';
        sampleFilename = 'P08samples.nii.gz';
    case "9"
        coreFilename = 'P09FLR.nii.gz';
        sampleFilename = 'P09T1G_points.nii.gz';
    case "10"
        coreFilename = 'P10FLR.nii.gz';
        sampleFilename = 'P10T1G_points.nii.gz';
    case "11"
        coreFilename = 'P11T1G.nii.gz';
        sampleFilename = 'P11T1G_points.nii.gz';
    case "12"
        coreFilename = 'P12FLR.nii.gz';
        sampleFilename = 'P12T1G_points.nii.gz';
    case "13"
        coreFilename = 'P13T1G.nii.gz';
        sampleFilename = 'P13T1G_points.nii.gz';
    case "14"
        coreFilename = 'P14FLR.nii.gz';
        sampleFilename = 'P14T1G_points.nii.gz';
    case "15"
        coreFilename = 'P15FLR.nii.gz';
        sampleFilename = 'P15T1G_points.nii.gz';
    case "17"
        coreFilename = 'P17T1G.nii.gz';
        sampleFilename = 'P17T1G_points.nii.gz';
end

% Plot tumor core volume
figure
core = extractCoords(coreFilename);
core.FaceColor = 'y';
core.FaceAlpha = 0.5;
core.EdgeColor = 'none';

% Plot corresponding tumor samples
hold on
[~,s] = extractSamples(sampleFilename);
samples = plot3(s(:,1),s(:,2),s(:,3),'ok');
samples.MarkerFaceColor = 'k'

% Calculate tumor distances, corresponding points on tumor surface,
% centroid, raidius, corresponding point on tumor surface, and volume
[distances, tumorpoints, C, R, Rpoint, vol] = calcDist(s,core);

% Plot tumor centroid and a line from the centroid to the furthest point on
% tumor surface (aka effective tumor radius)
centroid = plot3(C(1),C(2),C(3),'dr');
centroid.MarkerFaceColor = 'r';
centroid.MarkerSize = 10;

radius = plot3([C(1) Rpoint(1)],[C(2) Rpoint(2)],[C(3) Rpoint(3)],'-r');
radius.LineWidth = 2;

% Plot tumor distances with a symbol at the point on the tumor surfaces
% that the distance is calculated from
surfacepoints = plot3(tumorpoints(:,1),tumorpoints(:,2),tumorpoints(:,3),'bs');
surfacepoints.MarkerFaceColor = 'b';
surfacepoints.MarkerSize = 8;

distanceLine = plot3([s(:,1) tumorpoints(:,1)]',[s(:,2) tumorpoints(:,2)]',[s(:,3) tumorpoints(:,3)]','--b');

% Change Axis Font for clarity; add labels
set(gca,'linewidth',2)
set(gca,'fontsize',14)
set(get(gca,'XLabel'),'string','X (mm)')
set(get(gca,'YLabel'),'string','Y (mm)')
set(get(gca,'ZLabel'),'string','Z (mm)')

% Plot figure Legend, then adjust it to remove individual distance lines
legend('Tumor Core Volume', 'Biopsy Sample', 'Tumor Core Centroid', 'Effective Tumor Radius','Tumor Surface Point','Sample Distance')

end