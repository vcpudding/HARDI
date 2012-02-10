function compareDWIODFAnisotropicFiltering

close all
clear all

load GradientOrientations_64
UnitVectors
order = 4;
delta = 200;
G=constructMatrixOfMonomials(g,order);

SNRs = [10 25 50 100];
bVals = [1500 3000];
conductances = 1e-4;
weights = [0, 10, 100, 10000];

for iSNR=1:length(SNRs)
    for iBVal = 1:length(bVals)
        
        noisyImg = simulateCrossingFibers([16,16,3], 4, bVals(iBVal), SNRs(iSNR));
        noisyImg = permute(noisyImg, [4 1 2 3]);
        noisyODF = zeros(length(G), 16*16*3);

        for i=1:16*16*3            
            S = noisyImg(:,i);
            [ODF_coef, ODF_tensor] = Estimate_tensorODF(S, 1, GradientOrientations, order, delta);
            noisyODF(:,i) = G*ODF_coef;
            noisyODF(:,i) = noisyODF(:,i)/sum(noisyODF(:,i));
        end        
        
        for iCond = 1:length(conductances)
%             for iWeight = 1:length(weights)
%                 dwiFiltered = dwiAnisotropicFiltering(noisyImg, 1/SNRs(iSNR), conductances(iCond), weights(iWeight), 1e-3, 1e-6);                
%                 fileName = sprintf('results/DWI_filtered [SNR=%d][bVal=%d][c=%f][w=%d].mat', SNRs(iSNR), bVals(iBVal), conductances(iCond), weights(iWeight));
%                 save(fileName, 'dwiFiltered');
%                 disp(fileName);
%             end
            
            odfFiltered = odfAnisotropicFiltering(noisyODF, conductances(iCond), 1e-3, 1e-7);
            fileName = sprintf('results/ODF_filtered [SNR=%d][bVal=%d][c=%f].mat', SNRs(iSNR), bVals(iBVal), conductances(iCond));
            save(fileName, 'odfFiltered');
            disp(fileName);
        end
    end    
end

end