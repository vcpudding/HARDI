function filteredImg = anisotropicFiltering (noisyImg, conductance, tolerance)

filteredImg = noisyImg;
[m n] = size(noisyImg);
e = energyFunc(filteredImg, conductance);
eBuf = [];

while 1
    laste = e;
    step = 1e-5;
    dI = zeros(m,n);
    for i=1:m
        for j=1:n            
            gN = filteredImg(max(1,i-1),j) - filteredImg(i,j);
            gS = filteredImg(min(m,i+1),j) - filteredImg(i,j);
            gW = filteredImg(i, max(1,j-1)) - filteredImg(i,j);
            gE = filteredImg(i, min(n,j+1)) - filteredImg(i,j);
            dI(i,j) = exp(-gN^2/conductance)*gN+exp(-gS^2/conductance)*gS+exp(-gW^2/conductance)*gW+exp(-gE^2/conductance)*gE;
        end
    end
    
    filteredImg = filteredImg + step*dI;
    e = energyFunc(filteredImg, conductance);
    eBuf = [eBuf, e-laste];
    if abs(e-laste)<tolerance
        break;
    end
end

figure, plot(eBuf);

end

function e = energyFunc(img, conductance)

[gx gy] = gradient(img);
e = exp(-(gx.^2+gy.^2)/conductance);
e = conductance*sum(e(:));

end
