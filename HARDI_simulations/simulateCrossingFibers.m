function volDat = simulateCrossingFibers (dims, fiberWidth, bVal, snr)

load GradientOrientations_64
UnitVectors
order = 4;
delta = 200;

volDat = zeros([dims, length(GradientOrientations)]);

border1 = (dims(1)-fiberWidth)/2;
border2 = (dims(2)-fiberWidth)/2;

sig = 1/snr;

for i=1:dims(1)
    for j=1:dims(2)

        if i-border1>=1 && i-border1<=fiberWidth
            if j-border2>=1 && j-border2 <=fiberWidth
                %crossing region
                S = simulateDWData(bVal, GradientOrientations, [pi/2, pi], [0.5, 0.5], 1);
            else
                %horizontal region
                S = simulateDWData(bVal, GradientOrientations, 0, 1, 1);
            end
        else
            if j-border2>=1 && j-border2 <=fiberWidth
                %vertical region
                S = simulateDWData(bVal, GradientOrientations, pi/2, 1, 1);
            else
                %background
                S = simulateDWData(bVal, GradientOrientations, 0, 1, 0);
            end
        end
        
        for k=1:dims(3)
            volDat(i,j,k,:) = abs(S+sig*(randn(size(S))+sqrt(-1)*randn(size(S))));
        end
    end
end

% close all;
% v = reshape(volDat(:,:,2,:), [dims(1:2),length(GradientOrientations)]);
% sliceTitles = regexp(sprintf('(%0.2f, %0.2f, %0.2f)\n', GradientOrientations), '\n', 'split');
% plotVolumeData(v, 4, 4, sliceTitles);
% 
% odfField = zeros([dims, 15]);
% for i=1:dims(1)
%     for j=1:dims(2)
%         for k=1:dims(3)
%             S = reshape(volDat(i,j,1,:), [length(GradientOrientations),1]);
%             [ODF_coef, ODF_tensor] = Estimate_tensorODF(S, 1, GradientOrientations, order, delta);
%             odfField(i,j,k,:) = ODF_coef;
%         end
%     end
% end
% 
% plotODFField(odfField(:,:,2,:));
            

end