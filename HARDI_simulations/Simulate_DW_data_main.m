clear all
close all

fiber_direction1=[0 1 0]; % first fiber direction 
orientation=atan2(fiber_direction1(2),fiber_direction1(1));

angle=pi/4; % separation angle
b=3000; % acquisition b-value
w1=0.5; w2=0.5; % fiber weights
S0=1;
order=4;
delta=200;
sig=0.1; % SNR=1/sig
No_of_simulations=27;

R=[cos(-angle) sin(-angle) 0;-sin(-angle) cos(-angle) 0;0 0 1];
fiber_direction2=fiber_direction1*R'; % second fiber direction

load GradientOrientations_64
UnitVectors

% simulate data

S=Simulate_DW_data(b,GradientOrientations,orientation,angle,w1,w2);
%S = simulateDWData(b, GradientOrientations, [pi/2, 3*pi/4], [w1,w2], 1);
[ODF_coef ODF_tensor]=Estimate_tensorODF(S,S0,GradientOrientations,order,delta);
G=constructMatrixOfMonomials(g,order);
ODF = G*ODF_coef;
ODF = ODF/sum(ODF);
disp(GFA(ODF));

%figure,plotTensors(ODF_coef,0.9,[321 1 0]); % plot the mean ODF

for i=1:No_of_simulations

% corrupt by Rician noise
    y=randn(length(S),2);
    S_noisy_complex=S+sig*(y(:,1)+sqrt(-1)*y(:,2));
    S_noisy = abs(S_noisy_complex);

    %[ODF_coef ODF_tensor]=Estimate_tensorODF(S_noisy,S0,GradientOrientations,order,delta);
    
    %G=constructMatrixOfMonomials(g,order);
    %ODF=G*ODF_coef;
    
    %ODF=sqrt(ODF/sum(ODF(:))); % project ODF to the Hilbert sphere
    
    % compute ODF from coefficients
    
    DWI_matrix(:,i)=S_noisy;
    %sq_ODF_matrix(:,i)=ODF; 
    
end 

% mean_complex_DWI = mean(DWI_matrix,2);
% [ODF_coef ODF_tensor]=Estimate_tensorODF(abs(mean_complex_DWI),S0,GradientOrientations,order,delta);

% mean_ODF=calc_ODF_mean(sq_ODF_matrix);
% mean_ODF_coef=G\(mean_ODF.^2); % get the coefficients of the mean ODF
% figure,plotTensors(mean_ODF_coef,0.9,[321 1 0]); % plot the mean ODF

mean_DWI = removeRicianNoiseSingleVoxel(DWI_matrix, sig, 1e-8);
[ODF_coef ODF_tensor]=Estimate_tensorODF(mean_DWI,S0,GradientOrientations,order,delta);
figure,plotTensors(ODF_coef,0.9,[321 1 0]); % plot the mean ODF

ODF = G*ODF_coef;
ODF = ODF/sum(ODF);
disp(GFA(ODF));






    

