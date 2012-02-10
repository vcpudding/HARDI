clear noisyImg;
clear filteredImg;

noisyImg = simulateCrossingFibers([16,16,4], 4, 1500, 50);
filteredImg = dwiAnisotropicFiltering(noisyImg, 0.02, 1, 0, 1.0e-8);
gradIdx = 1;
figure, plotVolumeData(noisyImg(:,:,:,gradIdx), 1, 4);
figure, plotVolumeData(filteredImg(:,:,:, gradIdx), 1, 4);