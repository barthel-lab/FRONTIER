function [distances, tumorpoints, C, R, Rpoint, vol] = calcDist(sampleCoords,tumorvol)
%calcDist - calculates the shortest distance between the sample and surface
%of tumor volume
% sampleCoords - matrix of sample coordinates
% tumorvol - output patch object generated by extractCoords function

% Load vertices from patch object
v = tumorvol.Vertices;
f = tumorvol.Faces;

%%% Calculate centroid
% Rename corners
A = v(f(:,1),:);
B = v(f(:,2),:);
C = v(f(:,3),:);
% Needs to be **unnormalized** normals
N = cross(B-A,C-A,2);
% total volume via divergence theorem: ? 1
vol = sum(sum(v(f(:,1),:).*N))/6;
% centroid via divergence theorem and midpoint quadrature: ? x
C = 1/(2*vol)*(1/24* sum(N.*((A+B).^2 + (B+C).^2 + (C+A).^2)));

%%% Calculate the max Euclidean distance between centroid
%%% and a point on the surface of the tumor volume
for i = 1:size(v(:,1))
    if i == 1
        R = sqrt((v(i,1)-C(1))^2 + (v(i,2)-C(2))^2 + (v(i,3)-C(3))^2);
    else
        tempR = sqrt((v(i,1)-C(1))^2 + (v(i,2)-C(2))^2 + (v(i,3)-C(3))^2);
        if tempR > R
            R = tempR;
            Rpoint = v(i,:);
        end
    end
end

%%% Calculate the min Euclidean distance between sample coordinates
%%% and a point on the surface of the tumor volume
distances = zeros(size(sampleCoords(:,1))); % Array for euclidean distances
tumorpoints = zeros(size(sampleCoords)); % Array for the corresponding points on tumor volume surface

for i = 1:size(sampleCoords,1)
    for j = 1:size(v,1)
        if j == 1
            dist = sqrt((sampleCoords(i,1)-v(j,1))^2 + (sampleCoords(i,2)-v(j,2))^2 + (sampleCoords(i,3)-v(j,3))^2);
            tumorpoints(i,:) = v(j,:);
        else
            temp = sqrt((sampleCoords(i,1)-v(j,1))^2 + (sampleCoords(i,2)-v(j,2))^2 + (sampleCoords(i,3)-v(j,3))^2);
            if temp < dist
                dist = temp;
                tumorpoints(i,:) = v(j,:);
            end
        end
    end
    distances(i) = dist;
end

end