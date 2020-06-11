function status = extractCoords(filename)
%extractCoords - Extract coordinates from .nii file
%   View untouched .nii file using "Tools for NIfTI and ANALYZE image"
%   and extract coordinates from figure

P = load_untouch_nii(filename);

% Run package to produce interactive 3-View plot of tumor volume
%status = view_nii(file);  
%axesObjs = get(status.fig,'Children');
%dataObjs = get(axesObjs,'Children');

% Isolate image from loaded struct file
img = P.img;

% Plot isosurface of tumor volume
% img = smooth3(img); %  Smoothed isosurface
[f,v]=isosurface(img); % Extract patch faces and vertices data

%%%% Convert vox coordinates to real world space (mm) %%%%

% Load scaling settings from .nii header
% pixdim = P.hdr.dime.pixdim;
% qoffset_x = P.hdr.hist.qoffset_x;
% qoffset_y = P.hdr.hist.qoffset_y;
% qoffset_z = P.hdr.hist.qoffset_z;

srow_x = P.hdr.hist.srow_x;
srow_y = P.hdr.hist.srow_y;
srow_z = P.hdr.hist.srow_z;

%%% Matrices necessary to reorient and rescale voxels %%%

% qform reorientation matrix
% b = P.hdr.hist.quatern_b;
% c = P.hdr.hist.quatern_c;
% d = P.hdr.hist.quatern_d;
% a = sqrt(1-b^2-c^2-d^2);
% 
% R = [(a^2+b^2-c^2-d^2), 2*(b*c-a*d), 2*(b*d+a*c);
%      2*(b*c+a*d), (a^2-b^2+c^2-d^2), 2*(c*d-a*b); 
%      2*(b*d-a*c), 2*(c*d+a*b), (a^2-b^2-c^2+d^2)];
 
% sform affine matrix
S = [srow_x(1) srow_x(2) srow_x(3) srow_x(4);
     srow_y(1) srow_y(2) srow_y(3) srow_y(4);
     srow_z(1) srow_z(2) srow_z(3) srow_z(4);
         0         0         0         1    ];
     
%%% Transformation Calculations %%%

% Shorthand variables for cleaner equations
i = v(:,2);
j = v(:,1);
k = v(:,3);
%wq = zeros(size(v));
ws = zeros(size(v));

% Quaternion representation of rotation matrix (method 2 of NIfTI1 format)
% for ind = 1:length(i)
%     temp = R * [i(ind); j(ind); pixdim(1)*k(ind)] .* [pixdim(2); pixdim(3); pixdim(4)] + [qoffset_x; qoffset_y; qoffset_z];
%     wq(ind,1) = temp(1);
%     wq(ind,2) = temp(2);
%     wq(ind,3) = temp(3);
% end

% Affine transformation (method 3 of NIfTI1 format)
for ind = 1:length(i)
    temp = S * [i(ind); j(ind); k(ind); 1];
    ws(ind,1) = temp(1);
    ws(ind,2) = temp(2);
    ws(ind,3) = temp(3);
end

%%% Plot according to affine transformation method %%%
status = patch('Faces',f,'Vertices',ws,'FaceColor','yellow','FaceAlpha',0.3,'EdgeColor','none');


end

