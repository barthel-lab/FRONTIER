% Function to read the methylation class of a given datapoint and output
% the proper color for plotting
function class_color = findColor(input)
    switch input
        case 'Classic-like'
            class_color = 'r';
        case 'Codel'
            class_color = 'g';
        case 'Control'
            class_color = 'k';
        case 'G-CIMP-high'
            class_color = 'b';
        case 'G-CIMP-low'
            class_color = 'c';
        case 'LGm6-PA'
            class_color = 'm';
        case 'Mesenchymal-like'
            class_color = 'y';
    end
end