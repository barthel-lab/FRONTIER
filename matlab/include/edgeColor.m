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
