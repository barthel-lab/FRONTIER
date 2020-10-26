function wmat = plotWM()
hvoxcortfp = '/usr/local/fsl/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz'
hvoxsubfp = '/usr/local/fsl/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr0-1mm.nii.gz'

wmat = extractCoords(hvoxsubfp,[1,12])
wmat.FaceColor = '[0.25, 0.25, 0.25]';
wmat.EdgeColor = 'none';
wmat.FaceAlpha = 0.25;
end