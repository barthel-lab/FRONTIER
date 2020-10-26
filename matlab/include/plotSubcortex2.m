function [pall, hipp, amyg, accu] = plotSubcortex2()
hvoxcortfp = '/usr/local/fsl/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz'
hvoxsubfp = '/usr/local/fsl/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr0-1mm.nii.gz'

pall = extractCoords(hvoxsubfp,[7,18])
pall.FaceColor = '[0.8500, 0.3250, 0.0980]'; % orange / pallidum
pall.EdgeColor = 'none';
pall.FaceAlpha = 0.8;

hipp = extractCoords(hvoxsubfp,[9,19])
hipp.FaceColor = '[0.4940, 0.1840, 0.5560]'; % purple / hippocampus
hipp.EdgeColor = 'none';
hipp.FaceAlpha = 0.8;

amyg = extractCoords(hvoxsubfp,[10,20]) % blue / amygdala
amyg.FaceColor = '[0, 0, 1]';
amyg.EdgeColor = 'none';
amyg.FaceAlpha = 0.8

accu = extractCoords(hvoxsubfp,[11,21])
accu.FaceColor = '[1, 0, 0]'; % red / nucleus accumbens
accu.EdgeColor = 'none';
accu.FaceAlpha = 0.8;
end