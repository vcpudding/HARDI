function hFigs = plotODFField (odfCoefVol)

hFigs = [];
dims = size(odfCoefVol);
for k=1:dims(3)
    hFigs = [hFigs, figure('Position', [10,10,800,800])];
    for i=1:dims(1)
        for j=1:dims(2)
            subplot(dims(1), dims(2), j+dims(2)*(i-1));
            %subaxis(dims(1), dims(2), j+dims(2)*(i-1), 'Margin', 0);
            coef = reshape(odfCoefVol(i,j,k,:), [15,1]);
            plotTensors(coef, 0.9, [321, 0, -1]);
            axis off;
        end
    end
end

end