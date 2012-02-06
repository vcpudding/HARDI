function filteredImg = dwiAnisotropicFiltering (noisyImg, stdDeviation, conductance, priorWeight, tolerance)

%initialize
%filteredImg = vectorAvrFilter3D(noisyImg, 3);
filteredImg = noisyImg;
sqStdDeviation = stdDeviation^2;
logPost = logPosterior(filteredImg, noisyImg, stdDeviation, conductance, priorWeight);
a=0.001;
b=0.5;
it = 1;
eBuf = [];
%gradient descent
while 1
    lastLogPost = logPost;
    
    x = filteredImg.*noisyImg/sqStdDeviation;
    B = -noisyImg/sqStdDeviation + besseli(1,x,1)./besseli(0,x,1).*noisyImg/sqStdDeviation;
    
    P = vectorAnisotropicDivergence(filteredImg, conductance);
    
    step = 1e-3;
    %g = B+priorWeight*P;
    g = P;
    
    %disp(['step = ', num2str(step)]);
    filteredImg = filteredImg + step*g;
    
    logPost = logPosterior(filteredImg, noisyImg, stdDeviation, conductance, priorWeight);
    %disp(['posterior:', num2str(logPost)]);
    disp(['e:', num2str(logPost-lastLogPost)]);
    eBuf = [eBuf, logPost-lastLogPost];
    if abs(logPost-lastLogPost) <tolerance
        break;
    end;
    
    it = it+1;
    
end

%output
figure, plot(eBuf);

end

function p = logPosterior(filteredImg, noisyImg, sqStdDeviation, conductance, priorWeight)

z = filteredImg.*noisyImg/sqStdDeviation;
p1 = log(noisyImg/sqStdDeviation)-(filteredImg.^2+noisyImg.^2)/2/sqStdDeviation+log(besseli(0,z,1))+abs(z);
p1 = sum(p1(:));

[imgGradX, imgGradY, imgGradZ] = vectorGradient(filteredImg);
imgGradSqNorm = sum(imgGradX.^2,4)+sum(imgGradY.^2,4)+sum(imgGradZ.^2,4);
E = exp(-imgGradSqNorm/2/conductance);
p2 = priorWeight*conductance*sum(E(:));
%p2 = 1/(1-exp(-priorWeight*numel(noisyImg))) * exp(-E);

%p = p1+p2;
p=p2;

end

function [gx gy gz] = vectorGradient (vol)

vecdim = size(vol,4);
gx = zeros(size(vol));
gy = zeros(size(vol));
gz = zeros(size(vol));

for i=1:vecdim
    [g1 g2 g3] = gradient(vol(:,:,:,i));
    gx(:,:,:,i) = g1;
    gy(:,:,:,i) = g2;
    gz(:,:,:,i) = g3;
end

end


function div = vectorAnisotropicDivergence (vol, conductance)

[m n p vecdim] = size(vol);
div = zeros(size(vol));
for i=1:m
    for j=1:n
        for k=1:p
            gN = vol(max(1,i-1),j,k,:) - vol(i,j,k,:);
            gS = vol(min(m,i+1),j,k,:) - vol(i,j,k,:);
            gW = vol(i, max(1,j-1),k,:) - vol(i,j,k,:);
            gE = vol(i, min(n,j+1),k,:) - vol(i,j,k,:);
            gF = vol(i, j, max(1,k-1),:) - vol(i,j,k,:);
            gP = vol(i, j, min(p,k+1),:) - vol(i,j,k,:);
            
            div(i,j,k,:) = exp(-dot(gN,gN)/conductance)*gN + ...
                           exp(-dot(gS,gS)/conductance)*gS + ...
                           exp(-dot(gW,gW)/conductance)*gW + ...
                           exp(-dot(gE,gE)/conductance)*gE + ...
                           exp(-dot(gF,gF)/conductance)*gF + ...
                           exp(-dot(gP,gP)/conductance)*gP;
        end
    end
end
end