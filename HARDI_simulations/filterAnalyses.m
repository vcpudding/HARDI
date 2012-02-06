clc
clear all
close all

S0=1;
order=4;
delta=100;
load GradientOrientations_64

for iFile=1 %for all files
    
    fileName='';
    load(fileName); %load ODF_coef_mat, S_noisy_complex_mat and S_noisy
    
    %Filter complex signal
    S_complex = mean(S_noisy_complex, 1);
    [avrODF1 tensor1] = Estimate_tensorODF(abs(S_complex), S0, GradientOrientations, order, delta);
    
    %Filter magnitude signal
    S_magnitude % = S_noisy without Rician noise;
    [avrODF2 tensor2] = Estimate_tensorODF(S_magnitude, S0, GradientOrientations, order, delta);
    
    %Filter ODF
    ODFs %= Evaluate ODFs from tensor
    ODFs = sqrt(ODFs);%Take square roots
    avrODF3 = calc_ODF_mean(ODFs);
    
end

