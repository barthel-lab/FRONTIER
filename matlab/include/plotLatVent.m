function vent = plotLatVent()
    hvoxcortfp = '/usr/local/fsl/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz'
    hvoxsubfp = '/usr/local/fsl/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr0-1mm.nii.gz'

    vent = extractCoords(hvoxsubfp,[3,14])
    vent.FaceColor = '[0.9290, 0.6940, 0.1250]';
    vent.EdgeColor = 'none';
    vent.FaceAlpha = 0.8;
end