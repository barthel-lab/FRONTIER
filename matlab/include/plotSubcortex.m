function [thal, caud, puta] = plotSubcortex()
    hvoxcortfp = '/usr/local/fsl/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz'
    hvoxsubfp = '/usr/local/fsl/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr0-1mm.nii.gz'

    thal = extractCoords(hvoxsubfp,[4,15])
    thal.FaceColor = '[0, 0.5, 0]'; % green / thalamus
    thal.EdgeColor = 'none';
    thal.FaceAlpha = 0.8;

    caud = extractCoords(hvoxsubfp,[5,16])
    caud.FaceColor = '[0.3010, 0.7450, 0.9330]'; % light blue / caudatus
    caud.EdgeColor = 'none';
    caud.FaceAlpha = 0.8;

    puta = extractCoords(hvoxsubfp,[6,17])
    puta.FaceColor = '[1 0.4 0.6]'; % pink / putamen
    puta.EdgeColor = 'none';
    puta.FaceAlpha = 0.8;
end