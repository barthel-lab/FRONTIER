function distances = calcEucDist(point, points)
    %calcEucDist - calculates euclidian distances between a point and a set of
    %points
    
    % Instantiate output
    distances = zeros(size(points(:,1)));
    
    for i = 1:size(points,1)
        distances(i) = sqrt((points(i,1)-point(1))^2 + (points(i,2)-point(2))^2 + (points(i,3)-point(3))^2);
    end

end