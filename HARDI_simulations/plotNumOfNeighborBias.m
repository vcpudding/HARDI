close all
weight = [0.5,0.5];
dAngle = 45;
snr = 10;
bVal = 3000;

fileName = sprintf('SimData/Bias [Weight=(%0.1f, %0.1f)][sepAngle=%d][SNR=%d][bVal=%d].mat', weight(1), weight(2), dAngle, snr, bVal);
load(fileName);
sq_real_ODF = sqrt(real_ODF/sum(real_ODF));

nSamples = 100;
nNeighbors = 1:10:101;
biases = zeros(1, length(nNeighbors));

% load GradientOrientations_64
% UnitVectors
% order = 4;
% G=constructMatrixOfMonomials(g,order);

for i=1:length(nNeighbors)
    disp(['Neighbor no. = ', num2str(nNeighbors(i))]);
    odfLog = 0;
    for j=1:nSamples
        randIdx = randperm(size(sq_ODF_matrix,2));
        rand_sq_ODF_matrix = sq_ODF_matrix(:, randIdx(1:nNeighbors(i)));
        if nNeighbors(i)==1
            sq_mean_ODF = rand_sq_ODF_matrix(:,1);
        else
            sq_mean_ODF = calc_ODF_mean(rand_sq_ODF_matrix);
        end;
        %mean_ODF_coef = G\(mean_ODF.^2);
        %figure, plotTensors(mean_ODF_coef,0.9,[321 1 0]);
        d = dot(sq_mean_ODF, sq_real_ODF);
        l = (sq_real_ODF-d*sq_mean_ODF)/sqrt(1-d*d)*acos(d);
        odfLog = odfLog+l;
        if sum(isnan(odfLog))
            stop = 1;
        end
    end
    odfLog = odfLog/nSamples;
    disp(norm(odfLog));
    biases(i) = norm(odfLog);
end

figure;
h = plot(nNeighbors, biases);
saveas(h, 'comparison_output/neighbor_bias.fig');