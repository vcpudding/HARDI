function deviation = dirDeviation (odfCoef, realDirs)

r = size(realDirs,2);
odfTensor = convertCoef2Tensor(odfCoef);
P = cp_als(tensor(odfTensor), r, 'init', 'nvecs', 'tol', 1e-8, 'printitn', 0);
deviation = zeros(1,r);
for i=1:r
    tempDevs = acos(P{1}(:,i)'*realDirs);
    idx = tempDevs>pi/2;
    tempDevs(idx) = pi-tempDevs(idx);
    deviation(i) = min(tempDevs)*180/pi;
end

deviation = mean(deviation);

end