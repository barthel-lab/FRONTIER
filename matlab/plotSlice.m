function slice = plotSlice(img,plane,slicenum)
%plotSlice - isolates a 2D image from 3D .nii image matrix
%and plots that image in a 3D space
%   Function arguments:
%     img - the .nii image matrix
%     plane = axi/sag/cor - the plane that the 2D image will be isolated
%     from
%       axi - 
%       sag - the x-value corresponds to the slice number
%       cor - 
%     slicenum - integer for the desired slice (img size = 180x288x288)

% Use squeeze function to extract 2D slice from 3D image matrix
switch plane
    case "axi"
        slice = squeeze(img(:,:,slicenum));
    case "sag"
        origslice = squeeze(img(slicenum,:,:));
        % Reorient slice to match tumor space
        slice = imrotate(origslice,90);
    case "cor"
        slice = squeeze(img(:,slicenum,:));
end

% Lower image brightness threshold to show all of image
figure
imshow(slice,[0 500])
t = hgtransform('Parent',gca);
tx = makehgtform('translate',[0 slicenum 0]);
set(t,'Matrix',tx)

%ax = axes('XLim',[-37.1429 2.6182],'YLim',[16.4600 47.8200],'ZLim',[-26.9300 62.6700])

end

