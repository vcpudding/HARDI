function testODFAnisotropicFiltering (conductance, tolerance, step)

% noisyImg = simulateCrossingFibers([16,16,3], 4, 1500, 50);
% 
load GradientOrientations_64
UnitVectors
order = 4;
delta = 200;
G=constructMatrixOfMonomials(g,order);
% 
% noisyODF = zeros(16,16,3,length(G));
% 
% for i=1:16
%     for j=1:16
%         for k=1:3
%             S = reshape(noisyImg(i,j,1,:), [length(GradientOrientations),1]);
%             [ODF_coef, ODF_tensor] = Estimate_tensorODF(S, 1, GradientOrientations, order, delta);
%             noisyODF(i,j,k,:) = G*ODF_coef;
%             noisyODF(i,j,k,:) = noisyODF(i,j,k,:)/sum(noisyODF(i,j,k,:));
%         end
%     end
% end

load('../../Data/SimData/ODF_SNR=50.mat');
%noisyODF = permute(noisyODF, [4,1,2,3]);

filteredODF = odfAnisotropicFiltering(noisyODF, conductance, tolerance, step);
filteredODF = permute(filteredODF, [4 1 2 3]);
[n p] = size(filteredODF);
filteredODF_coef = G\reshape(filteredODF, [n p]);
filteredODF_coef = ipermute(filteredODF_coef, [4 1 2 3]);
figure, plotODFField(filteredODF_coef(:,:,2,:));

end