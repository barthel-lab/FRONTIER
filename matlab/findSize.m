

function marker_size = findSize(input)
    switch string(input)
        case 'TRUE'
            marker_size = 10;
        case 'FALSE'
            marker_size = 5;
    end
end