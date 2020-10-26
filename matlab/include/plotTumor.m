function plotTumor(patient)
    % Nii file map
    niimapfp = '../results/meta/FRONTIER.nii.transformed.csv'

    % Nii file directory
    filedir = '../data/masks/nii2_mni';
    
    % Load nii filemap
    niifiles = readtable(niimapfp)
    sniifiles = niifiles(find(strcmp(niifiles.Patient,patient)),:)
    
    % Prepare nii file paths
    t1gfile = ''
    if(strlength(sniifiles.t1g) > 0)
        t1gfile = char(strcat(filedir, '/', sniifiles.t1g));
    end
    flairfile = ''
    if(strlength(sniifiles.flair) > 0)
        flairfile = char(strcat(filedir, '/', sniifiles.flair))
    end
    maskfile = char(strcat(filedir, '/', sniifiles.mask))
    
    % Draw T1G volume
    if(strlength(t1gfile) > 0)
        tumor1 = extractCoords(t1gfile, 65000, false);
        tumor1.FaceColor = 'red'
        tumor1.EdgeColor = 'none';
        tumor1.FaceAlpha = 0.20;
    end
    
    % Draw FLAIR volume
    if(strlength(flairfile) > 0)
        tumor2 = extractCoords(flairfile, 65000, false);
        tumor2.FaceColor = [1 0.65 0];
        tumor2.EdgeColor = 'none';
        tumor2.FaceAlpha = 0.15;
    end
end