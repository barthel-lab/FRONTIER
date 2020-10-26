function cort = plotCortex()
    hvoxcortfp = '/usr/local/fsl/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz'
    hvoxsubfp = '/usr/local/fsl/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr0-1mm.nii.gz'

    cort = extractCoords(hvoxcortfp)
    cort.FaceColor = 'black';
    cort.EdgeColor = 'none';
    cort.FaceAlpha = 0.08;
end