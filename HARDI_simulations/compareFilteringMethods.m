function compareFilteringMethods(weight, snr, bVal)

close all;
addpath('tensor_toolbox_2.4');
addpath('tensor_toolbox_2.4/algorithms');
sepAngles = [45, 60, 75, 90];

load GradientOrientations_64
UnitVectors
order=4;
delta=100;
G=constructMatrixOfMonomials(g,order);
S0 = 1;

tableDistFR = cell(4,5);
tableDistFR(1,:) = {'Fisher Rao distance', '45', '60', '75', '90'};
tableDistFR(2:4, 1) = {'Mean complex DWIs','Mean magnitude DWIs','FR mean of ODFs'};
tableDirDev = cell(5, 5);
tableDirDev(1,:) = {'Direction deviation', '45', '60', '75', '90'};
tableDirDev(2:5, 1) = {'Mean complex DWIs','Mean magnitude DWIs','FR mean of ODFs', 'Real DWI'};

for iAngle = 1:length(sepAngles)
    
    fileName = sprintf('SimData/Data [Weight=(%0.1f, %0.1f)][sepAngle=%d][SNR=%d][bVal=%d].mat', weight(1), weight(2), sepAngles(iAngle), snr, bVal);
    load(fileName);    
    
    fiber_direction1=[0 1 0]; % first fiber direction 
    angle = sepAngles(iAngle)*pi/180;
    orientation=atan2(fiber_direction1(2),fiber_direction1(1));
    R=[cos(-angle) sin(-angle) 0;-sin(-angle) cos(-angle) 0;0 0 1];
    fiber_direction2=fiber_direction1*R'; % second fiber direction
    fiberDirs = [fiber_direction1', fiber_direction2'];

    real_ODF=G*real_ODF_coef; 
    complexDWI_ODF=G*mean_ODF_coef_avrComplexDWI; 
    magnitudeDWI_ODF=G*mean_ODF_coef_avrMagnitudeDWI; 
    FR_ODF=G*mean_ODF_coef_FR; 
    
    tableDistFR{2,iAngle+1} = fisherRaoDist(complexDWI_ODF, real_ODF);
    tableDistFR{3,iAngle+1} = fisherRaoDist(magnitudeDWI_ODF, real_ODF);
    tableDistFR{4,iAngle+1} = fisherRaoDist(FR_ODF, real_ODF);
    
    tableDirDev{2,iAngle+1} = dirDeviation(mean_ODF_coef_avrComplexDWI, fiberDirs);
    tableDirDev{3,iAngle+1} = dirDeviation(mean_ODF_coef_avrMagnitudeDWI, fiberDirs);
    tableDirDev{4,iAngle+1} = dirDeviation(mean_ODF_coef_FR, fiberDirs);
    tableDirDev{5,iAngle+1} = dirDeviation(real_ODF_coef, fiberDirs);
    
end

disp(tableDistFR);
disp(tableDirDev);


end