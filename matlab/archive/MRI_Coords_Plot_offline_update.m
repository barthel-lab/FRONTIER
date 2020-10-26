% MRI_Coords_Plot_offline.m - plots the coordinates for MRI biopsy data
% Script assumes metadata excel sheet has already been loaded to workspace
% as a table named VUmc

if ~exist('VUmc')
    tablefile = '../sandbox/Sample_Sheet_VUmc_Location_Stats.xls';
    VUmc = readtable(tablefile);
end

setup = input('Load patient data? (Y/N) ','s');
% Only run if data isn't loaded into workspace
if strcmpi(setup,'y') || strcmpi(setup,'yes')
    % Identify VUmc Patient Samples (designated with a number from 1 to 17)
    patients = unique(VUmc.Patient);

    % Assumes VUmc patients column is imported as double class
    p1 = find(strcmp(VUmc.Patient,'VUmc-01'));
    p2 = find(strcmp(VUmc.Patient,'VUmc-02'));
    p3 = find(strcmp(VUmc.Patient,'VUmc-03'));
    p4 = find(strcmp(VUmc.Patient,'VUmc-04'));
    p5 = find(strcmp(VUmc.Patient,'VUmc-05'));
    p6 = find(strcmp(VUmc.Patient,'VUmc-06'));
    p7 = find(strcmp(VUmc.Patient,'VUmc-07'));
    p8 = find(strcmp(VUmc.Patient,'VUmc-08'));
    p9 = find(strcmp(VUmc.Patient,'VUmc-09'));
    p10 = find(strcmp(VUmc.Patient,'VUmc-10'));
    p11 = find(strcmp(VUmc.Patient,'VUmc-11'));
    p12 = find(strcmp(VUmc.Patient,'VUmc-12'));
    p13 = find(strcmp(VUmc.Patient,'VUmc-13'));
    p14 = find(strcmp(VUmc.Patient,'VUmc-14'));
    p15 = find(strcmp(VUmc.Patient,'VUmc-15'));
    p17 = find(strcmp(VUmc.Patient,'VUmc-17'));

    % Create subset tables for each patient
    Tp1 = VUmc(p1,:);
    Tp2 = VUmc(p2,:);
    Tp3 = VUmc(p3,:);
    Tp4 = VUmc(p4,:);
    Tp5 = VUmc(p5,:);
    Tp6 = VUmc(p6,:);
    Tp7 = VUmc(p7,:);
    Tp8 = VUmc(p8,:);
    Tp9 = VUmc(p9,:);
    Tp10 = VUmc(p10,:);
    Tp11 = VUmc(p11,:);
    Tp12 = VUmc(p12,:);
    Tp13 = VUmc(p13,:);
    Tp14 = VUmc(p14,:);
    Tp15 = VUmc(p15,:);
    Tp17 = VUmc(p17,:);
end


p = input('Enter desired patient dataset (format: Tp[x]): ');
figure
hold on
for i = 1:size(p,1)
    % switch loop identifies the location specification for a given
    % datapoint and plots the corresponding symbol. The appropriate colors
    % are chosen by the findColor, edgeColor functions
    switch p.Location{i,1}
        case 'C'
            class_color = findColor(p.Cell_Predict{i,1}); % Determine marker color
            s = strcat('o',class_color); % Combine marker color with symbol spec
            tx_color = edgeColor(p.Tx2017_Predict{i,1}); % Determine marker edge color
            % For Batch 1, used X/Y/Z_corrected; for Batch 2, use X/Y/Z
            %if p.Batch(1) == '1'
                plot3(p.X_corrected(i),p.Y_corrected(i),p.Z_corrected(i),s,'MarkerFaceColor',class_color,'MarkerEdgeColor',tx_color,'MarkerSize',15,'LineWidth',2.5)
            %else
            %    plot3(p.X(i),p.Y(i),p.Z(i),s,'MarkerFaceColor',class_color,'MarkerEdgeColor',tx_color,'MarkerSize',15,'LineWidth',2.5)
            %end
        case 'P'
            class_color = findColor(p.Cell_Predict{i,1});
            s = strcat('h',class_color);
            tx_color = edgeColor(p.Tx2017_Predict{i,1}); 
            %if p.Batch(1) == '1'
                plot3(p.X_corrected(i),p.Y_corrected(i),p.Z_corrected(i),s,'MarkerFaceColor',class_color,'MarkerEdgeColor',tx_color,'MarkerSize',15,'LineWidth',2.5)
            %else
            %    plot3(p.X(i),p.Y(i),p.Z(i),s,'MarkerFaceColor',class_color,'MarkerEdgeColor',tx_color,'MarkerSize',15,'LineWidth',2.5)
            %end
        case 'No'
            class_color = findColor(p.Cell_Predict{i,1});
            s = strcat('^',class_color);
            tx_color = edgeColor(p.Tx2017_Predict{i,1}); 
            %if p.Batch(1) == '1'
                plot3(p.X_corrected(i),p.Y_corrected(i),p.Z_corrected(i),s,'MarkerFaceColor',class_color,'MarkerEdgeColor',tx_color,'MarkerSize',15,'LineWidth',2.5)
            %else
            %    plot3(p.X(i),p.Y(i),p.Z(i),s,'MarkerFaceColor',class_color,'MarkerEdgeColor',tx_color,'MarkerSize',15,'LineWidth',2.5)
            %end
        case 'NA'
            class_color = findColor(p.Cell_Predict{i,1});
            s = strcat('*',class_color);
            tx_color = edgeColor(p.Tx2017_Predict{i,1}); 
            %if p.Batch(1) == '1'
                plot3(p.X_corrected(i),p.Y_corrected(i),p.Z_corrected(i),s,'MarkerFaceColor',class_color,'MarkerEdgeColor',tx_color,'MarkerSize',15,'LineWidth',2.5)
            %else
            %    plot3(p.X(i),p.Y(i),p.Z(i),s,'MarkerFaceColor',class_color,'MarkerEdgeColor',tx_color,'MarkerSize',15,'LineWidth',2.5)
            %end
    end
    if i == size(p,1)
        axis equal
    end
end


% Change Axis Font for clarity; add labels
set(gca,'linewidth',2)
set(gca,'fontsize',14)
set(get(gca,'XLabel'),'string','X (mm)')
set(get(gca,'YLabel'),'string','Y (mm)')
set(get(gca,'ZLabel'),'string','Z (mm)')

% Nii file directory
filedir = '../data/nii';

% Plot tumor volume(s) on same figure
% T1G volumes are yellow, FLR volumes are magenta'
switch p.Patient{1}
    case "VUmc-01"
        filename = ['P01T1G.nii.gz';'P01FLR.nii.gz'];
        color = ['m';'y'];
    case "VUmc-02"
        filename = ['P02T1G.nii.gz';'P02FLR.nii.gz'];
        color = ['m';'y'];
    case "VUmc-03"
        filename = 'P03_new_volume.nii.gz'; % T1G
        color = 'm';
    case "VUmc-04"
        filename = ['P04T1G.nii.gz';'P04FLR.nii.gz'];
        color = ['m';'y'];
    case "VUmc-05"
        filename = 'P05T1G.nii.gz';
        color = 'm';
    case "VUmc-06"
        filename = 'P06T1G.nii.gz';
        color = 'm';
    case "VUmc-07"
        filename = ['P07T1G.nii.gz';'P07FLR.nii.gz'];
        color = ['m';'y'];
    case "VUmc-08"
        filename = ['P08T1G.nii.gz';'P08FLR.nii.gz'];
        color = ['m';'y'];
    case "VUmc-09"
        filename = 'P09FLR.nii.gz';
        color = 'y';
    case "VUmc-10"
        filename = 'P10FLR.nii.gz';
        color = 'y';
    case "VUmc-11"
        filename = ['P11T1G.nii.gz';'P11FLR.nii.gz'];
        color = ['m';'y'];
    case "VUmc-12"
        filename = 'P12FLR.nii.gz';
        color = 'y';
    case "VUmc-13"
        filename = ['P13T1G.nii.gz';'P13FLR.nii.gz'];
        color = ['m';'y'];
    case "VUmc-14"
        filename = 'P14FLR.nii.gz';
        color = 'y';
    case "VUmc-15"
        filename = 'P15FLR.nii.gz';
        color = 'y';
    case "VUmc-17"
        filename = ['P17T1G.nii.gz';'P17FLR.nii.gz'];
        color = ['m';'y'];
end

if size(filename,1) == 2 
    tumor1 = extractCoords(strcat(filedir,'/',filename(1,:)));
    tumor1.FaceColor = color(1);
    tumor1.EdgeColor = 'none';
    tumor1.FaceAlpha = 0.20;
    tumor2 = extractCoords(strcat(filedir,'/',filename(2,:)));
    tumor2.FaceColor = color(2);
    tumor2.EdgeColor = 'none';
    tumor2.FaceAlpha = 0.20;
    patient = ['Patient',' ',p.Patient{1}];
    set(get(gca,'title'),'string',patient)
else
    tumor = extractCoords(strcat(filedir,'/',filename));
    tumor.FaceColor = color;
    tumor.EdgeColor = 'none';
    tumor.FaceAlpha = 0.20;
    patient = ['Patient',' ',p.Patient{1}];
    set(get(gca,'title'),'string',patient)
end

%% Add brain overlay
hold on
brain_overlay = extractCoords('../data/nii/P17_brain.mask.nii.gz')
brain_overlay.FaceColor = 'black';
brain_overlay.EdgeColor = 'none';
brain_overlay.FaceAlpha = 0.08;


% Render mp4 video of rotating 3D figure, if desired
video = input('Create rotating video? (Y/N) ','s');
if strcmpi(video,'y') || strcmpi(video,'yes')
    filename = input('Enter filename: ','s');
    OptionZ.FrameRate=5;OptionZ.Duration=15;OptionZ.Periodic=true; % Rendering settings
    CaptureFigVid([-20,10;-110,35;-190,80;-290,10;-380,10],filename,OptionZ)
end

% Plot legend on separate figure using placeholder points
legendPlot = input('Plot legend? (Y/N) ','s'); % Run only if desired
if strcmpi(legendPlot,'y') || strcmpi(legendPlot,'yes')
    figure
    hold on
    % Location Symbols
    plot(NaN,NaN,'.w','MarkerFaceColor','w','MarkerSize',12) % Subheader
    plot(NaN,NaN,'ok','MarkerFaceColor','k','MarkerSize',12) % Core
    plot(NaN,NaN,'hk','MarkerFaceColor','k','MarkerSize',12) % Periphery
    plot(NaN,NaN,'^k','MarkerFaceColor','k','MarkerSize',12) % Outside
    plot(NaN,NaN,'*k','MarkerFaceColor','k','MarkerSize',12) % N/A
    plot(NaN,NaN,'.w','MarkerFaceColor','w','MarkerSize',12) % Filler Space
    % Methylation Classification Colors
    plot(NaN,NaN,'.w','MarkerFaceColor','w','MarkerSize',12) % Subheader
    plot(NaN,NaN,'sr','MarkerFaceColor','r','MarkerSize',12) % Classic-like
    plot(NaN,NaN,'sg','MarkerFaceColor','g','MarkerSize',12) % Codel
    plot(NaN,NaN,'sk','MarkerFaceColor','k','MarkerSize',12) % Cortex
    plot(NaN,NaN,'sb','MarkerFaceColor','b','MarkerSize',12) % G-CIMP-high
    plot(NaN,NaN,'sc','MarkerFaceColor','c','MarkerSize',12) % Inflammatory-TME
    plot(NaN,NaN,'sm','MarkerFaceColor','m','MarkerSize',12) % Reactive-TME
    plot(NaN,NaN,'sy','MarkerFaceColor','y','MarkerSize',12) % Mesenchymal-like
    plot(NaN,NaN,'.w','MarkerFaceColor','w','MarkerSize',12) % Filler Space
    % Transcription Classification Colors
    plot(NaN,NaN,'.w','MarkerFaceColor','w','MarkerSize',12) % Subheader
    plot(NaN,NaN,'-','Color',[0.85 0.33 0.1],'LineWidth',2.0) % Classical
    plot(NaN,NaN,'-','Color',[0.47 0.67 0.19],'LineWidth',2.0) % Mesenchymal
    plot(NaN,NaN,'-','Color',[0.49 0.18 0.56],'LineWidth',2.0) % Proneural
    % Legend Labels
    legend('Location Symbol:','Core','Periphery','Outside Tumor','Unspecified Location','','Classification Color:','Classic-like','Codel','Control','G-CIMP-high','G-CIMP-low','LGm6-PA','Mesenchymal-like','','Transcription Subtype Color:','Classical','Mesenchymal','Proneural')
end

%%%%% Subfunctions %%%%%

% Function to read the methylation class of a given datapoint and output
% the proper color for plotting
function class_color = findColor(input)
    switch input
        case 'Classic-like'
            class_color = 'r';
        case 'Codel'
            class_color = 'g';
        case 'Cortex'
            class_color = 'k';
        case 'G-CIMP-high'
            class_color = 'b';
        case 'Inflammatory-TME'
            class_color = 'c';
        case 'Reactive-TME'
            class_color = 'm';
        case 'Mesenchymal-like'
            class_color = 'y';
    end
end

% Function to read the transcription class of a given datapoint and output
% the proper color for plotting
function tx_color = edgeColor(input)
    switch input
        case 'CL'
            tx_color = [0.85 0.33 0.1];
        case 'MS'
            tx_color = [0.47 0.67 0.19];
        case 'PN'
            tx_color = [0.49 0.18 0.56];
    end
end

