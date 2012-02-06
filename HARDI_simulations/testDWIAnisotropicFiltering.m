clear noisyImg;
clear filteredImg;

noisyImg = simulateCrossingFibers([16,16,4], 4, 1500, 50);

gradIdx = 1;
filteredImg = dwiAnisotropicFiltering(noisyImg, 0.1, 100, 1, 1.0e-8);
figure(h1);
plotVolumeData(noisyImg(:,:,:,gradIdx), 1, 4);
figure(h2);
plotVolumeData(filteredImg(:,:,:, gradIdx), 1, 4);