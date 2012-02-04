clc
clear all
close all

weights = [0.3, 0.4, 0.5];
sepAngles = [45, 60, 75, 90];
SNRs = [5,10,25,50,100];
bVals = [1500, 3000];

load GradientOrientations_64
UnitVectors
order=4;
delta=100;
G=constructMatrixOfMonomials(g,order);

for iWeights = 1:length(weights)
    for iAngles = 1:length(sepAngles)
        for iSNRs = 1:length(SNRs)
            for iBVals = 1:length(bVals)                

                dAngle = sepAngles(iAngles);
                angle=dAngle*pi/180; % separation angle
                b=bVals(iBVals); % acquisition b-value
                w1=weights(iWeights); w2=1-weights(iWeights); % fiber weights
                S0=1;
                snr = SNRs(iSNRs);
                sig=1/snr; % SNR=1/sig

                fiber_direction1=[0 1 0]; % first fiber direction 
                orientation=atan2(fiber_direction1(2),fiber_direction1(1));
                R=[cos(-angle) sin(-angle) 0;-sin(-angle) cos(-angle) 0;0 0 1];
                fiber_direction2=fiber_direction1*R'; % second fiber direction

                % simulate data
                S=Simulate_DW_data(b,GradientOrientations,orientation,angle,w1,w2);
                
                for iNoise = 1:100
                    % corupt by Rician noise
                    y=randn(length(S),2);
                    S_noisy_complex = S+sig*(y(:,1)+sqrt(-1)*y(:,2));
                    S_noisy = abs(S_noisy_complex);
                    
                    [ODF_coef tensor]=Estimate_tensorODF(S_noisy,S0,GradientOrientations,order,delta);
    
                    ODF=G*ODF_coef;
                    ODF=sqrt(ODF/sum(ODF(:))); % project ODF to the Hilbert sphere

                    % compute ODF from coefficients

                    DWI_matrix(:,iNoise)=S_noisy;
                    complex_DWI_matrix(:,iNoise) = S_noisy_complex;
                    sq_ODF_matrix(:,iNoise)=ODF; 

                end  
                
                [real_ODF_coef real_tensor] = Estimate_tensorODF(S,S0,GradientOrientations,order,delta);
                
                mean_DWI_complex = mean(complex_DWI_matrix, 2);
                [mean_ODF_coef_avrComplexDWI mean_tensor_avrComplexDWI] = Estimate_tensorODF(abs(mean_DWI_complex),S0,GradientOrientations,order,delta);
                
                mean_DWI_magnitude = removeRicianNoiseSingleVoxel(DWI_matrix,sig, 1e-8);
                [mean_ODF_coef_avrMagnitudeDWI mean_tensor_avrMagnitudeDWI] = Estimate_tensorODF(mean_DWI_magnitude,S0,GradientOrientations,order,delta);
                
                mean_ODF_FR=calc_ODF_mean(sq_ODF_matrix);
                mean_ODF_coef_FR=G\(mean_ODF_FR.^2); % get the coefficients of the mean ODF     

                fileName = sprintf('SimData/Data [Weight=(%0.1f, %0.1f)][sepAngle=%d][SNR=%d][bVal=%d].mat', w1, w2, dAngle, snr, b);
                display(fileName);
                save(fileName, 'S', 'DWI_matrix', 'complex_DWI_matrix', 'sq_ODF_matrix', 'real_ODF_coef', 'mean_ODF_coef_avrComplexDWI', 'mean_ODF_coef_avrMagnitudeDWI', 'mean_ODF_coef_FR');
                
            end
        end
    end    
end