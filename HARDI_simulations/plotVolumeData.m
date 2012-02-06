function hFigs = plotVolumeData (volDat, rowsPerFig, colsPerFig, sliceTitles)

nSlicePerFig = rowsPerFig*colsPerFig;
dims = size(volDat);

hFigs = [];
for i=1:dims(3)
    if mod(i, nSlicePerFig)==1
        %hFigs = [hFigs, figure];
        colormap gray;
        ratio = length(colormap)/max(max(max(volDat)));
    end
    subplot(rowsPerFig, colsPerFig, mod(i-1, nSlicePerFig)+1);    
    borderedSlice = zeros(dims(1:2)+2);
    borderedSlice(2:dims(1)+1, 2:dims(2)+1) = volDat(:,:,i)*ratio;
    image(borderedSlice);
    if nargin>3
        title(sliceTitles{i});
    end
    axis square;
    axis off;
end

end