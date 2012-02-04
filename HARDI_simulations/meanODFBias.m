function meanODFBias (weight, snr, bVal)

sepAngles = [45, 60, 75, 90];

load GradientOrientations_64
UnitVectors
order=4;
delta=100;
G=constructMatrixOfMonomials(g,order);
nSimulations = 100000;

for iAngle = 1:length(sepAngles)
    disp (['angle = ', num2str(sepAngles(iAngle))]);
    
    dAngle = sepAngles(iAngle);
    angle=dAngle*pi/180; % separation angle
    b=bVal; % acquisition b-value
    w1=weight(1); w2=weight(2); % fiber weights
    S0=1;
    sig=1/snr; % SNR=1/sig

    fiber_direction1=[0 1 0]; % first fiber direction 
    orientation=atan2(fiber_direction1(2),fiber_direction1(1));
    R=[cos(-angle) sin(-angle) 0;-sin(-angle) cos(-angle) 0;0 0 1];
    fiber_direction2=fiber_direction1*R'; % second fiber direction

    % simulate data
    S=Simulate_DW_data(b,GradientOrientations,orientation,angle,w1,w2);

    [real_ODF_coef real_tensor] = Estimate_tensorODF(S,S0,GradientOrientations,order,delta);
    real_ODF = G*real_ODF_coef;
    real_ODF = real_ODF/sum(real_ODF);
    
    ODF_coef_matrix = zeros(size(G,2), nSimulations);
    noiseReal = randn(length(S), nSimulations);
    noiseImag = randn(length(S), nSimulations);
    S_noisy_complex = S*ones(1, nSimulations)+sig*(noiseReal+sqrt(-1)*noiseImag);
    S_noisy = abs(S_noisy_complex);
    for iNoise = 1:nSimulations
        [ODF_coef ODF_tensor]=Estimate_tensorODF(S_noisy(:, iNoise),S0,GradientOrientations,order,delta);        
        ODF_coef_matrix(:, iNoise) = ODF_coef;
        if mod(iNoise,100)==0
            disp(['sim #', num2str(iNoise)]);
        end
    end

    ODF_matrix = G*ODF_coef_matrix;
    sum_ODF_matrix = ones(size(ODF_matrix,1), 1)*sum(ODF_matrix, 1);
    sq_ODF_matrix = sqrt(ODF_matrix./sum_ODF_matrix);
        
    fileName = sprintf('SimData/Bias [Weight=(%0.1f, %0.1f)][sepAngle=%d][SNR=%d][bVal=%d].mat', w1, w2, dAngle, snr, b);
    display(fileName);
    save(fileName, 'nSamples', 'real_ODF', 'sq_ODF_matrix');
end

end