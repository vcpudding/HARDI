function plotAnisotropicFilteringResults

close all
clear all

addpath('tensor_toolbox_2.4');
addpath('tensor_toolbox_2.4/algorithms');

SNRs = [10 25 50 100];
bVals = [1500 3000];
conductances = [1e-4, 1e-3, 1e-2];
weights = [0, 10, 100, 10000];
dims = [16,16,3];
fiberWidth = 4;

for iBVal = 1%:length(bVals)        
    %cleanImg = simulateCrossingFibers(dims, fiberWidth, bVals(iBVal), 10000);
    %cleanODFSlice = odfFromImg(cleanImg(:,:,2,:));
    load('testCleanODFSlice.mat');
        
    for iSNR=1%:length(SNRs)
        
        for iWeight = 1%:length(weights)
            dFRs = zeros(dims(1), dims(2), length(conductances));
            dirDevs = zeros(dims(1), dims(2), length(conductances));
            gfas = zeros(dims(1), dims(2), length(conductances));
            for iCond = 1:length(conductances)
                fileName = sprintf('results/DWI_filtered [SNR=%d][bVal=%d][c=%f][w=%d].mat', SNRs(iSNR), bVals(iBVal), conductances(iCond), weights(iWeight));
                load(fileName);
                dwiFilteredODFSlice = odfFromImg(dwiFiltered(:,:,2,:));
                [dFR dirDev gfa] = stats(dwiFilteredODFSlice, cleanODFSlice, fiberWidth);
                dFRs(:,:,iCond) = dFR;
                dirDevs(:,:,iCond) = dirDev;
                gfas(:,:,iCond) = gfa;
                if iCond>1
                    continue;
                end
                plotAndSaveODF(dwiFilteredODFSlice, sprintf('figures/DWI_odf|SNR=%d|bVal=%d|c=%g|w=%d.eps', SNRs(iSNR), bVals(iBVal), conductances(iCond), weights(iWeight)));
            end
        
            plotAndSaveFigure(dFRs, [min(dFRs(:)), max(dFRs(:))], sprintf('figures/DWI_dFR|SNR=%d|bVal=%d|w=%d.eps', SNRs(iSNR), bVals(iBVal), weights(iWeight)));
            plotAndSaveFigure(dirDevs, [0, 20], sprintf('figures/DWI_dirDev|SNR=%d|bVal=%d|w=%d.eps', SNRs(iSNR), bVals(iBVal), weights(iWeight)));
            plotAndSaveFigure(gfas, [min(gfas(:)), max(gfas(:))], sprintf('figures/DWI_GFA|SNR=%d|bVal=%d|w=%d.eps', SNRs(iSNR), bVals(iBVal), weights(iWeight)));
        end
        
        dFRs = zeros(dims(1), dims(2), length(conductances));
        dirDevs = zeros(dims(1), dims(2), length(conductances));
        gfas = zeros(dims(1), dims(2), length(conductances));
        for iCond = 1:length(conductances)        
            fileName = sprintf('results/ODF_filtered [SNR=%d][bVal=%d][c=%f].mat', SNRs(iSNR), bVals(iBVal), conductances(iCond));
            load(fileName);
            
            [dFR dirDev gfa] = stats(odfFiltered(:,:,2,:), cleanODFSlice, fiberWidth);
            dFRs(:,:,iCond) = dFR;
            dirDevs(:,:,iCond) = dirDev;
            gfas(:,:,iCond) = gfa;
            if iCond>1
                continue;
            end
            plotAndSaveODF(odfFiltered(:,:,2,:), sprintf('figures/ODF_odf|SNR=%d|bVal=%d|c=%g.eps', SNRs(iSNR), bVals(iBVal), conductances(iCond)));
        end
        
        plotAndSaveFigure(dFRs, [min(dFRs(:)), max(dFRs(:))], sprintf('figures/ODF_dFR|SNR=%d|bVal=%d.eps', SNRs(iSNR), bVals(iBVal)));
        plotAndSaveFigure(dirDevs, [0, 20], sprintf('figures/ODF_dirDev|SNR=%d|bVal=%d.eps', SNRs(iSNR), bVals(iBVal)));
        plotAndSaveFigure(gfas, [min(gfas(:)), max(gfas(:))], sprintf('figures/ODF_GFA|SNR=%d|bVal=%d.eps', SNRs(iSNR), bVals(iBVal)));
    end    
end


end

function [dFR dirDev gfa] = stats (odfSlice, cleanODFSlice, fiberWidth)

[m n odfDim] = size(odfSlice);

dFR = zeros(m,n);
dirDev = zeros(m,n);
gfa = zeros(m,n);

load GradientOrientations_64
UnitVectors
order = 4;
G=constructMatrixOfMonomials(g,order);

for i=1:m
    for j=1:n
        odf = reshape(odfSlice(i,j,:), [odfDim,1]);
        cleanODF = reshape(cleanODFSlice(i,j,:), [odfDim, 1]);
        dFR(i,j) = fisherRaoDist(odf, cleanODF);
        odfCoef = G\odf;
        dirDev(i,j) = dirDeviation(odfCoef, getRealFibDirs([m,n], fiberWidth,i,j));
        gfa(i,j) = GFA(odf)-GFA(cleanODF);
    end
end

end

function dirs = getRealFibDirs (dims, fiberWidth, i,j)
border1 = (dims(1)-fiberWidth)/2;
border2 = (dims(2)-fiberWidth)/2;

if i-border1>=1 && i-border1<=fiberWidth
    if j-border2>=1 && j-border2 <=fiberWidth
        %crossing region
        dirs = [0 1 0; 1 0 0]';
    else
        %horizontal region
        dirs = [1 0 0]';
    end
else
    if j-border2>=1 && j-border2 <=fiberWidth
        %vertical region
        dirs = [0 1 0]';
    else
        %background
        dirs = [1 0 0]';
    end
end

end

function plotAndSaveFigure (dat, colAxis, fileName)

hFig = figure('Visible', 'off','Position', [10 10 800 200]);
colormap jet;
caxis(colAxis);
nSubplots = size(dat,3);
ax = zeros(nSubplots,1);
for i=1:nSubplots
    ax(i) = subaxis(1,nSubplots,i,'Spacing',0,'ML',0,'Padding',0);
    image(dat(:,:,i), 'CDataMapping', 'scaled');
    axis square; axis off;
end
hBar=colorbar;
set(hBar, 'Position', [.9 0.34 .03 .32]);

for i=1:nSubplots
    pos=get(ax(i), 'Position');
    set(ax(i), 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
end
print(hFig, fileName, '-depsc');
end

function plotAndSaveODF (odf, fileName)

load GradientOrientations_64
UnitVectors
order = 4;
G=constructMatrixOfMonomials(g,order);

dims = size(odf);
odf = permute(odf, [4 1 2 3]);
odf = reshape(odf, [dims(4), prod(dims(1:3))]);
odfCoef = G\odf;
odfCoef = reshape(odfCoef, [size(odfCoef,1), dims(1:3)]);
odfCoef = ipermute(odfCoef, [4 1 2 3]);
hFig = plotODFField(odfCoef);
print(hFig, fileName, '-depsc');

end