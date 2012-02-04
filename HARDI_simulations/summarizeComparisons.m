clc
clear all
close all

weights = [0.3, 0.4, 0.5];
sepAngles = [45, 60, 75, 90];
SNRs = [5,10,25,50,100];
bVals = [1500, 3000];

addpath('tensor_toolbox_2.4');
addpath('tensor_toolbox_2.4/algorithms');

load GradientOrientations_64
UnitVectors
order=4;
delta=100;
G=constructMatrixOfMonomials(g,order);

copyfile('comparison_output/template.tex', 'comparison_output/comparisons.tex');
fid = fopen('comparison_output/comparisons.tex', 'a');
fprintf(fid, '\\begin{document}\n');

fprintf(fid, '{\\color{red} Red: Mean complex DWI}\\\\');
fprintf(fid, '{\\color{black} Black: Estimated magnitude DWI}\\\\');
fprintf(fid, '{\\color{blue} Blue: Mean ODF}\\\\');

%weight-angle
for iSNR = 1:length(SNRs)
    for iBVal = 1:length(bVals)
        snr = SNRs(iSNR);
        bVal = bVals(iBVal);

        for iWeights = 1:length(weights)
            for iAngle = 1:length(sepAngles)    
                weight = weights(iWeights);
                angle = sepAngles(iAngle);
                fileName = sprintf('SimData/Data [Weight=(%0.1f, %0.1f)][sepAngle=%d][SNR=%d][bVal=%d].mat', weight, 1-weight, angle, snr, bVal);
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

                frDistBuf(iWeights, iAngle, :) = [fisherRaoDist(complexDWI_ODF, real_ODF), fisherRaoDist(magnitudeDWI_ODF, real_ODF),fisherRaoDist(FR_ODF, real_ODF)];
                dirDevBuf(iWeights, iAngle, :) = [dirDeviation(mean_ODF_coef_avrComplexDWI, fiberDirs), dirDeviation(mean_ODF_coef_avrMagnitudeDWI, fiberDirs), dirDeviation(mean_ODF_coef_FR, fiberDirs)];
            end
        end

        fprintf(fid, '\\begin{table}[H]\n\\caption{SNR=%d, $b$-value=%d}\n\\begin{center}\n', snr, bVal);

        fprintf(fid, '\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}c |*{4}{c}}\n');
        fprintf(fid, '\\multicolumn{5}{c}{\\textbf{Fisher-Rao distances to the real ODF ($\\times 10^{-2}$)}}\\\\ \\hline\n');
        fprintf(fid, '\\backslashbox{Weights}{Separating angles} & $45^{\\circ}$ & $60^{\\circ}$ & $75^{\\circ}$ & $90^{\\circ}$ \\\\ \\hline\n');
        for iWeight = 1:length(weights)
            fprintf(fid, '(%0.1f, %0.1f)', weights(iWeight), 1-weights(iWeight));
            for iAngle = 1:length(sepAngles)
                 fprintf(fid, '& {\\color{red} %0.2f}\\;\\;{\\color{black} %0.2f}\\;\\;{\\color{blue} %0.2f}', frDistBuf(iWeight, iAngle, :)*1e2);
            end
            fprintf(fid, '\\\\\n');
        end
        fprintf(fid, '\\hline\n\\end{tabular*}\n');

        fprintf(fid, '\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}c |*{4}{c}}\n');
        fprintf(fid, '\\multicolumn{5}{c}{\\textbf{Direction deviation}}\\\\ \\hline\n');
        fprintf(fid, '\\backslashbox{Weights}{Separating angles} & $45^{\\circ}$ & $60^{\\circ}$ & $75^{\\circ}$ & $90^{\\circ}$ \\\\ \\hline\n');
        for iWeight = 1:length(weights)
            fprintf(fid, '(%0.1f, %0.1f)', weights(iWeight), 1-weights(iWeight));
            for iAngle = 1:length(sepAngles)
                 fprintf(fid, '& {\\color{red} %0.2f}\\;\\;{\\color{black} %0.2f}\\;\\;{\\color{blue} %0.2f}', dirDevBuf(iWeight, iAngle, :));
            end
            fprintf(fid, '\\\\\n');
        end
        fprintf(fid, '\\hline\n\\end{tabular*}\n');

        fprintf(fid, '\\end{center}\n\\end{table}\n\n\n');
        
    end
    
    fprintf(fid, '\\clearpage');
end

fprintf(fid, '\\end{document}');
fclose(fid);

return;

