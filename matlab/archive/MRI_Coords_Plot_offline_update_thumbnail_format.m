% MRI_Coords_Plot_offline.m - plots the coordinates for MRI biopsy data
% Script assumes metadata excel sheet has already been loaded to workspace
% as a table named VUmc

if ~exist('VUmc')
    tablefile = '/Users/anderk/Documents/FRONTIER/Metadata/Sample_Sheet_VUmc_Location_Stats-082819.txt';
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
% Plot sample biopsies with sample number
plot3(p.X_corrected,p.Y_corrected,p.Z_corrected,'ko','MarkerFaceColor','k','MarkerEdgeColor','none','MarkerSize',1,'LineWidth',0.25)
text(p.X_corrected,p.Y_corrected,p.Z_corrected,[repmat(' ',[length(p.Sample_No) 1]) num2str(p.Sample_No)],'HorizontalAlignment','left','FontSize',5,'FontName','Helvetica');



% Remove axis ticks and rescale figure
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
set(gca,'box','off','XColor','none','YColor','none','ZColor','none','TickDir','out')
set(gcf,'Position',[0 0 50 50])
set(gca,'Visible','off')

% Nii file directory
filedir = '/Users/anderk/Documents/FRONTIER/VUmc Patient MRI Coords/Nii_Files';

% Plot tumor volume(s) on same figure
% T1G volumes are yellow, FLR volumes are magenta'
switch p.Patient{1}
    case "VUmc-01"
        filename = ['P01T1G.nii.gz';'P01FLR.nii.gz'];
        color = [[1 0 1];[1 0.84 0]]; % color = [magenta; gold]
        az = -69; % Azimuth angle of view
        el = 6; % Elevation angle of view
        % Orientation arrow vectors origin (subscript "o") and endpoint
        % (subscript "e")
        xo = [-50 60 50];
        xe = [-30 60 50];
        yo = [-50 60 50];
        ye = [-50 80 50];
        zo = [-50 60 50];
        ze = [-50 60 70];
    case "VUmc-02"
        filename = ['P02T1G.nii.gz';'P02FLR.nii.gz'];
        color = [[1 0 1];[1 0.84 0]];
        az = 94;
        el = 11;
        xo = [40 -30 20];
        xe = [60 -30 20];
        yo = [40 -30 20];
        ye = [40 -10 20];
        zo = [40 -30 20];
        ze = [40 -30 40];
    case "VUmc-03"
        filename = 'P03_new_FLR.nii.gz';
        color = [1 0.84 0];
        az = -158;
        el = 14;
        xo = [40 50 50];
        xe = [50 50 50];
        yo = [40 50 50];
        ye = [40 60 50];
        zo = [40 50 50];
        ze = [40 50 60];
    case "VUmc-04"
        filename = ['P04T1G.nii.gz';'P04FLR.nii.gz'];
        color = [[1 0 1];[1 0.84 0]];
        az = -138;
        el = -2;
        xo = [-18 -5 5];
        xe = [-8 -5 5];
        yo = [-18 -5 5];
        ye = [-18 5 5];
        zo = [-18 -5 5];
        ze = [-18 -5 15];
    case "VUmc-05"
        filename = 'P05FLR.nii.gz';
        color = [1 0.84 0];
        az = 115;
        el = 20;
        xo = [55 -25 5];
        xe = [65 -25 5];
        yo = [55 -25 5];
        ye = [55 -15 5];
        zo = [55 -25 5];
        ze = [55 -25 15];
    case "VUmc-06"
        filename = 'P06FLR.nii.gz';
        color = [1 0.84 0];
        az = -26;
        el = 35;
        xo = [-10 -10 -50];
        xe = [0 -10 -50];
        yo = [-10 -10 -50];
        ye = [-10 0 -50];
        zo = [-10 -10 -50];
        ze = [-10 -10 -40];
    case "VUmc-07"
        filename = ['P07T1G.nii.gz';'P07FLR.nii.gz'];
        color = [[1 0 1];[1 0.84 0]];
        az = -169;
        el = 11;
        xo = [-15 -35 -10];
        xe = [-5 -35 -10];
        yo = [-15 -35 -10];
        ye = [-15 -25 -10];
        zo = [-15 -35 -10];
        ze = [-15 -35 0];
    case "VUmc-08"
        filename = ['P08T1G.nii.gz';'P08FLR.nii.gz'];
        color = [[1 0 1];[1 0.84 0]];
        az = 28;
        el = 5;
        xo = [55 10 50];
        xe = [65 10 50];
        yo = [55 10 50];
        ye = [55 20 50];
        zo = [55 10 50];
        ze = [55 10 60];
    case "VUmc-09"
        filename = 'P09FLR.nii.gz';
        color = [1 0.84 0];
        az = 106;
        el = 26;
        xo = [-10 -15 60];
        xe = [-5 -15 60];
        yo = [-10 -15 60];
        ye = [-10 -10 60];
        zo = [-10 -15 60];
        ze = [-10 -15 65];
    case "VUmc-10"
        filename = 'P10FLR.nii.gz';
        color = [1 0.84 0];
        az = 151;
        el = 14;
        xo = [30 50 -20];
        xe = [40 50 -20];
        yo = [30 50 -20];
        ye = [30 60 -20];
        zo = [30 50 -20];
        ze = [30 50 -10];
    case "VUmc-11"
        filename = ['P11T1G.nii.gz';'P11FLR.nii.gz'];
        color = [[1 0 1];[1 0.84 0]];
        az = -138;
        el = 37;
        xo = [5 -55 0];
        xe = [15 -55 0];
        yo = [5 -55 0];
        ye = [5 -45 0];
        zo = [5 -55 0];
        ze = [5 -55 10];
    case "VUmc-12"
        filename = 'P12FLR.nii.gz';
        color = [1 0.84 0];
        az = -172;
        el = 15;
        xo = [-55 40 -5];
        xe = [-45 40 -5];
        yo = [-55 40 -5];
        ye = [-55 50 -5];
        zo = [-55 40 -5];
        ze = [-55 40 5];
    case "VUmc-13"
        filename = ['P13T1G.nii.gz';'P13FLR.nii.gz'];
        color = [[1 0 1];[1 0.84 0]];
        az = -120;
        el = 8;
        xo = [15 -50 -10];
        xe = [25 -50 -10];
        yo = [15 -50 -10];
        ye = [15 -40 -10];
        zo = [15 -50 -10];
        ze = [15 -50 0];
    case "VUmc-14"
        filename = 'P14FLR.nii.gz';
        color = [1 0.84 0];
        az = -60;
        el = 34;
        xo = [-35 -50 15];
        xe = [-25 -50 15];
        yo = [-35 -50 15];
        ye = [-35 -40 15];
        zo = [-35 -50 15];
        ze = [-35 -50 25];
    case "VUmc-15"
        filename = 'P15FLR.nii.gz';
        color = [1 0.84 0];
        az = -4;
        el = 20;
        xo = [-10 -5 -30];
        xe = [0 -5 -30];
        yo = [-10 -5 -30];
        ye = [-10 5 -30];
        zo = [-10 -5 -30];
        ze = [-10 -5 -20];
    case "VUmc-17"
        filename = ['P17T1G.nii.gz';'P17FLR.nii.gz'];
        color = [[1 0 1];[1 0.84 0]];
        az = -47;
        el = 35;
        xo = [-5 -35 -50];
        xe = [5 -35 -50];
        yo = [-5 -35 -50];
        ye = [-5 -25 -50];
        zo = [-5 -35 -50];
        ze = [-5 -35 -40];
end

if size(filename,1) == 2 
    tumor1 = extractCoords(strcat(filedir,'/',filename(1,:)));
    tumor1.FaceColor = color(1,:);
    tumor1.EdgeColor = 'none';
    tumor1.FaceAlpha = 0.10;
    tumor2 = extractCoords(strcat(filedir,'/',filename(2,:)));
    tumor2.FaceColor = color(2,:);
    tumor2.EdgeColor = 'none';
    tumor2.FaceAlpha = 0.10;
%     patient = ['Patient',' ',p.Patient{1}];
%     set(get(gca,'title'),'string',patient)
    view(az,el)
    axis equal
    hold on
    xarrow = mArrow3(xo,xe,'color','red');
    yarrow = mArrow3(yo,ye,'color','blue');
    zarrow = mArrow3(zo,ze,'color','green');
    hold off
else
    tumor = extractCoords(strcat(filedir,'/',filename));
    tumor.FaceColor = color;
    tumor.EdgeColor = 'none';
    tumor.FaceAlpha = 0.10;
%     patient = ['Patient',' ',p.Patient{1}];
%     set(get(gca,'title'),'string',patient)
    view(az,el)
    axis equal
    hold on
    xarrow = mArrow3(xo,xe,'color','red');
    yarrow = mArrow3(yo,ye,'color','blue');
    zarrow = mArrow3(zo,ze,'color','green');
    hold off
end


% Save figure to file
% set(gcf, 'PaperUnits', 'normalized')
% set(gcf, 'PaperPosition', [0 0 1 1])
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'-painters','-dpdf','-r300',['/Users/anderk/Documents/FRONTIER/VUmc Patient MRI Coords/Matlab_Figures/Thumbnails/',p.Patient{1},'.pdf'])

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

