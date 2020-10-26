% Function to read the methylation class of a given datapoint and output
% the proper color for plotting
function class_color = findColor(input)
    switch string(input)
        case 'Classic-like'
            class_color = [1 0 0];
        case 'Codel'
            class_color = [0.63 0.13 0.94];
        case 'Inflammatory-TME'
            class_color = [1 0.65 0];
        case 'G-CIMP-high'
            class_color = [0 1 0];
        case 'G-CIMP-low'
            class_color = [0 1 1];
        case 'LGm6-PA'
            class_color = [1 0 1];
        case 'Mesenchymal-like'
            class_color = [0 0 1];
        case 'Cortex'
            class_color = [0.55 0.35 0.17];
        case 'Reactive-TME'
            class_color = [1 1 0];
        case 'NA'
            class_color = [0 0 0];
    end
    class_color;
end