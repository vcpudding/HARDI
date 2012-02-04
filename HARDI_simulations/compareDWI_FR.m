function compareDWI_FR (weight, sepAngle, snr, bVal)

%load GradientOrientations_64
UnitVectors
order=4;
G=constructMatrixOfMonomials(g,order);

fileName = sprintf('SimData/Data [Weight=(%0.1f, %0.1f)][sepAngle=%d][SNR=%d][bVal=%d].mat', weight(1), weight(2), sepAngle, snr, bVal);
load(fileName);

real_ODF=G*real_ODF_coef; real_ODF = real_ODF/sum(real_ODF);
mean_ODF_avrDWI=G*mean_ODF_coef_avrComplexDWI; mean_ODF_avrDWI = mean_ODF_avrDWI/sum(mean_ODF_avrDWI);
mean_ODF_FR=G*mean_ODF_coef_FR; mean_ODF_FR = mean_ODF_FR/sum(mean_ODF_FR);

GFA_real = GFA(real_ODF);
GFA_avrDWI = GFA(mean_ODF_avrDWI);
GFA_FR = GFA(mean_ODF_FR);

dFR_avrDWI = acos(sqrt(mean_ODF_avrDWI)'*sqrt(real_ODF));
dFR_FR = acos(sqrt(mean_ODF_FR)'*sqrt(real_ODF));

figure;
subplot(1,3,1); plotTensors(real_ODF_coef, 0.9, [321,1,0]);title(sprintf('Real ODF\nGFA=%f', GFA_real)); 
subplot(1,3,2); plotTensors(mean_ODF_coef_avrComplexDWI, 0.9, [321,1,0]);title(sprintf('Gaussian mean of noisy DWIs\nGFA=%f\nd_{FR}=%f', GFA_avrDWI, dFR_avrDWI)); 
subplot(1,3,3); plotTensors(mean_ODF_coef_FR, 0.9, [321,1,0]);title(sprintf('FR mean of noisy ODFs\nGFA=%f\nd_{FR}=%f', GFA_FR, dFR_FR));

end