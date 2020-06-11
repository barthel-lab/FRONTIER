function [P,convertedCoords] = extractSamples(filename)
%extractSamples - extracts sample coordinates from MRI .nii file

% Load unscaled and unrotated .nii image
P = load_untouch_nii(filename);
img = P.img;

% Determine the indices corresponding to samples, which have max intensity
ind = find(img == 65535);
[i,j,k] = ind2sub(size(img),ind);

% Plot indices to see coordinates in voxel space
% plot3(i,j,k,'.');

% For a given sample, plot shows a diamond organization of 12 points. The
% sample coordinate is the centroid of those points. Z is constant for a
% given set of points, and given the point spacing, the centroid is simply
% [mean(x-axis points), mean(y-axis points), z]
%%% NOTE: This method does not work for patients 7,8; the 12 points are not
%%% grouped together for 2 sample coordinates
samplenum = length(i)/12;
coords = zeros(samplenum,3);
for x = 1:samplenum
    startind = 1+(12*(x-1));
    endind = x*12;
    coords(x,1) = mean(i(startind:endind));
    coords(x,2) = mean(j(startind:endind));
    coords(x,3) = k(startind);
end

% Initialization for coordinate space conversion
convertedCoords = zeros(samplenum,3);
i = coords(:,1);
j = coords(:,2);
k = coords(:,3);

% Perform affine transformation on coordinates to convert to real world
% space

% Load scaling settings from .nii header
srow_x = P.hdr.hist.srow_x;
srow_y = P.hdr.hist.srow_y;
srow_z = P.hdr.hist.srow_z;

% sform affine matrix
S = [srow_x(1) srow_x(2) srow_x(3) srow_x(4);
     srow_y(1) srow_y(2) srow_y(3) srow_y(4);
     srow_z(1) srow_z(2) srow_z(3) srow_z(4);
         0         0         0         1    ];

% Affine transformation (method 3 of NIfTI1 format)
for ind = 1:length(coords)
    temp = S * [i(ind); j(ind); k(ind); 1];
    convertedCoords(ind,1) = temp(1);
    convertedCoords(ind,2) = temp(2);
    convertedCoords(ind,3) = temp(3);
end

end