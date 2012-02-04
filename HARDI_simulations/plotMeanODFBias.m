function plotMeanODFBias (weight, snr, bVal, dAngle)

close all
fileName = sprintf('SimData/Bias [Weight=(%0.1f, %0.1f)][sepAngle=%d][SNR=%d][bVal=%d].mat', weight(1), weight(2), dAngle, snr, bVal);
load(fileName);
sq_real_ODF = sqrt(real_ODF/sum(real_ODF));

nSamples = 100:50:3700;
nNeighbors = 27;
biases = zeros(1, length(nSamples));

% for i=1:length(nSamples)
%     disp(['Sample no. = ', num2str(nSamples(i))]);
%     odfLog = 0;
%     for j=1:nSamples(i)
%         randIdx = randperm(size(sq_ODF_matrix,2));
%         rand_sq_ODF_matrix = sq_ODF_matrix(:, randIdx(1:nNeighbors));
%         sq_mean_ODF = calc_ODF_mean(rand_sq_ODF_matrix);
%         %mean_ODF_coef = G\(mean_ODF.^2);
%         %figure, plotTensors(mean_ODF_coef,0.9,[321 1 0]);
%         d = dot(sq_mean_ODF, sq_real_ODF);
%         l = (sq_real_ODF-d*sq_mean_ODF)/sqrt(1-d*d)*acos(d);
%         odfLog = odfLog+l;
%         %disp(norm(l));
%     end
%     odfLog = odfLog/nSamples(i);
%     biases(i) = norm(odfLog);
% end
% 
% figure;
% h1 = plot(nSamples, biases);
% saveas(h1, 'comparison_output/mean_ODF_bias.fig');

nSamples = 3000;
nNeighbors = 1:100;
biases = zeros(1, length(nNeighbors));
gfas = zeros(1, length(nNeighbors));

for i=1:length(nNeighbors)
    disp(['Neighbor no. = ', num2str(nNeighbors(i))]);
    odfLog = 0;
    gfa = 0;
    for j=1:nSamples
        randIdx = randperm(size(sq_ODF_matrix,2));
        rand_sq_ODF_matrix = sq_ODF_matrix(:, randIdx(1:nNeighbors(i)));
        if nNeighbors(i)==1
            sq_mean_ODF = rand_sq_ODF_matrix(:,1);
        else
            sq_mean_ODF = calc_ODF_mean(rand_sq_ODF_matrix);
        end;
        d = dot(sq_mean_ODF, sq_real_ODF);
        l = (sq_real_ODF-d*sq_mean_ODF)/sqrt(1-d*d)*acos(d);
        odfLog = odfLog+l;
        gfa = gfa+GFA(sq_mean_ODF.^2);
    end
    odfLog = odfLog/nSamples;
    gfa = gfa/nSamples;
    disp(norm(odfLog));
    biases(i) = norm(odfLog);
    gfas(i) = gfa;
end

figure;
h2 = plot(nNeighbors, biases);
saveas(h2, 'comparison_output/neighbor_bias.fig');

figure;
h3 = plot(nNeighbors, gfas);
saveas(h3, 'comparison_output/mean_ODF_gfa.fig');

end