% MRI_Coords_Plot_plotly.m - plots the coordinates for MRI biopsy data
% Script assumes metadata excel sheet has already been loaded to workspace
% as a table named MSG

setup = input('Load patient data? (Y/N) ','s');
% Only run if data isn't loaded into workspace
if strcmpi(setup,'y') || strcmpi(setup,'yes')
    % Identify VUmc Patient Samples (designated with a number from 1 to 8)
    patients = unique(MSG.Patient);

    % Assumes VUmc patients column is imported as double class
    p1 = find(MSG.Patient == 1);
    p2 = find(MSG.Patient == 2);
    p3 = find(MSG.Patient == 3);
    p4 = find(MSG.Patient == 4);
    p5 = find(MSG.Patient == 5);
    p6 = find(MSG.Patient == 6);
    p7 = find(MSG.Patient == 7);
    p8 = find(MSG.Patient == 8);
    p9 = find(MSG.Patient == 9);
    p10 = find(MSG.Patient == 10);
    p11 = find(MSG.Patient == 11);
    p12 = find(MSG.Patient == 12);
    p13 = find(MSG.Patient == 13);
    p14 = find(MSG.Patient == 14);
    p15 = find(MSG.Patient == 15);
    p17 = find(MSG.Patient == 17);

    % Create subset tables for each patient
    Tp1 = MSG(p1,:);
    Tp2 = MSG(p2,:);
    Tp3 = MSG(p3,:);
    Tp4 = MSG(p4,:);
    Tp5 = MSG(p5,:);
    Tp6 = MSG(p6,:);
    Tp7 = MSG(p7,:);
    Tp8 = MSG(p8,:);
    Tp9 = MSG(p9,:);
    Tp10 = MSG(p10,:);
    Tp11 = MSG(p11,:);
    Tp12 = MSG(p12,:);
    Tp13 = MSG(p13,:);
    Tp14 = MSG(p14,:);
    Tp15 = MSG(p15,:);
    Tp17 = MSG(p17,:);
end

p = input('Enter desired patient dataset (format: Tp[x]): ');
figure
hold on
for i = 1:size(p,1)
    % switch loop identifies the location specification for a given
    % datapoint and plots the corresponding symbol. The appropriate colors
    % are chosen by the findColor, edgeColor functions
    switch p.Location(i,1)
        case 'C'
            class_color = findColor(p.Cell_Predict(i,1)); % Determine marker color
            s = strcat('o',class_color); % Combine marker color with symbol spec
            tx_color = edgeColor(p.Tx_Predict(i,1)); % Determine marker edge color
            % For Batch 1, used X/Y/Z_corrected; for Batch 2, use X/Y/Z
            if p.Batch(1) == '1'
                scatter3(p.X_match(i),p.Y_match(i),p.Z_match(i),10,'MarkerFaceColor',class_color,'MarkerEdgeColor',tx_color,'linewidth',6)
            else
                scatter3(p.X(i),p.Y(i),p.Z(i),10,'MarkerFaceColor',class_color,'MarkerEdgeColor',tx_color,'linewidth',6)
            end
        case 'P'
            class_color = findColor(p.Cell_Predict(i,1));
            s = strcat('d',class_color);
            tx_color = edgeColor(p.Tx_Predict(i,1)); 
            if p.Batch(1) == '1'
                scatter3(p.X_match(i),p.Y_match(i),p.Z_match(i),10,'MarkerFaceColor',class_color,'MarkerEdgeColor',tx_color,'linewidth',6)
            else
                scatter3(p.X(i),p.Y(i),p.Z(i),10,'MarkerFaceColor',class_color,'MarkerEdgeColor',tx_color,'linewidth',6)
            end
        case 'No'
            class_color = findColor(p.Cell_Predict(i,1));
            s = strcat('s',class_color);
            tx_color = edgeColor(p.Tx_Predict(i,1)); 
            if p.Batch(1) == '1'
                scatter3(p.X_match(i),p.Y_match(i),p.Z_match(i),10,'MarkerFaceColor',class_color,'MarkerEdgeColor',tx_color,'linewidth',6)
            else
                scatter3(p.X(i),p.Y(i),p.Z(i),10,'MarkerFaceColor',class_color,'MarkerEdgeColor',tx_color,'linewidth',6)
            end
        case 'NA'
            class_color = findColor(p.Cell_Predict(i,1));
            s = strcat('*',class_color);
            tx_color = edgeColor(p.Tx_Predict(i,1)); 
            if p.Batch(1) == '1'
                scatter3(p.X_match(i),p.Y_match(i),p.Z_match(i),10,'MarkerFaceColor',class_color,'MarkerEdgeColor',tx_color,'linewidth',6)
            else
                scatter3(p.X(i),p.Y(i),p.Z(i),10,'MarkerFaceColor',class_color,'MarkerEdgeColor',tx_color,'linewidth',6)
            end
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

% Plot tumor and brain volume on same figure
switch p.Patient(1)
    case "1"
        filename = 'P01T1G.nii.gz';
    case "2"
        filename = 'P02T1G.nii.gz';
    case "3"
        filename = 'P03_new_volume.nii.gz';
    case "4"
        filename = 'P04T1G.nii.gz';
    case "5"
        filename = 'P05T1G.nii.gz';
    case "6"
        filename = 'P06T1G.nii.gz';
    case "7"
        filename = 'P07T1G.nii.gz';
    case "8"
        filename = 'P08T1G.nii.gz';
    case "9"
        filename = 'P09FLR.nii.gz';
    case "10"
        filename = 'P10FLR.nii.gz';
    case "11"
        filename = 'P11T1G.nii.gz';
    case "12"
        filename = 'P12FLR.nii.gz';
    case "13"
        filename = 'P13T1G.nii.gz';
    case "14"
        filename = 'P14FLR.nii.gz';
    case "15"
        filename = 'P15FLR.nii.gz';
    case "17"
        filename = 'P17T1G.nii.gz';
end

tumor = extractCoords(filename);
tumor.FaceColor = 'red';
tumor.EdgeColor = 'red';
tumor.FaceAlpha = 0.20;
patient = ['Patient',' ',num2str(p.Patient(1))];
set(get(gca,'title'),'string',patient)

% Render plotly interactive figure, if desired
interaction = input('Create plotly figure? (Y/N) ','s');
if strcmpi(interaction,'y') || strcmpi(interaction,'yes')
    response = fig2plotly(gcf,'strip',false);
end

% Plot legend on separate figure using placeholder points
legendPlot = input('Plot legend? (Y/N) ','s'); % Run only if desired
if strcmpi(legendPlot,'y') || strcmpi(legendPlot,'yes')
    figure
    hold on
    % Location Symbols
    plot(NaN,NaN,'.w','MarkerFaceColor','w','MarkerSize',12) % Subheader
    plot(NaN,NaN,'ok','MarkerFaceColor','k','MarkerSize',12) % Core
    plot(NaN,NaN,'dk','MarkerFaceColor','k','MarkerSize',12) % Periphery
    plot(NaN,NaN,'sk','MarkerFaceColor','k','MarkerSize',12) % Outside
    plot(NaN,NaN,'*k','MarkerFaceColor','k','MarkerSize',12) % N/A
    plot(NaN,NaN,'.w','MarkerFaceColor','w','MarkerSize',12) % Filler Space
    % Methylation Classification Colors
    plot(NaN,NaN,'.w','MarkerFaceColor','w','MarkerSize',12) % Subheader
    plot(NaN,NaN,'sr','MarkerFaceColor','r','MarkerSize',12) % Classic-like
    plot(NaN,NaN,'sg','MarkerFaceColor','g','MarkerSize',12) % Codel
    plot(NaN,NaN,'sk','MarkerFaceColor','k','MarkerSize',12) % Control
    plot(NaN,NaN,'sb','MarkerFaceColor','b','MarkerSize',12) % G-CIMP-high
    plot(NaN,NaN,'sc','MarkerFaceColor','c','MarkerSize',12) % G-CIMP-low
    plot(NaN,NaN,'sm','MarkerFaceColor','m','MarkerSize',12) % LGm6-PA
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
            class_color = 'red';
        case 'Codel'
            class_color = 'green';
        case 'Control'
            class_color = 'black';
        case 'G-CIMP-high'
            class_color = 'blue';
        case 'G-CIMP-low'
            class_color = 'cyan';
        case 'LGm6-PA'
            class_color = 'magenta';
        case 'Mesenchymal-like'
            class_color = 'yellow';
    end
end

% Function to read the transcription class of a given datapoint and output
% the proper color for plotting
function tx_color = edgeColor(input)
    switch input
        case 'CL'
            tx_color = uint8([0.85 0.33 0.1].*255 + 0.5);
        case 'MS'
            tx_color = uint8([0.47 0.67 0.19].*255 + 0.5);
        case 'PN'
            tx_color = uint8([0.49 0.18 0.56].*255 + 0.5);
    end
end

