function filteredImg = dwiAnisotropicFiltering (noisyImg, stdDeviation, conductance, priorWeight, tolerance)

%initialize
filteredImg = vectorAvrFilter3D(noisyImg, 3);
sqStdDeviation = stdDeviation^2;
logPost = logPosterior(filteredImg, noisyImg, stdDeviation, conductance, priorWeight);
a=0.1;
b=0.5;
it = 1;

%gradient descent
while 1
    lastLogPost = logPost;
    
    x = filteredImg.*noisyImg/sqStdDeviation;
    B = -noisyImg/sqStdDeviation + besseli(1,x,1)./besseli(0,x,1).*noisyImg/sqStdDeviation;
    [imgGradX, imgGradY, imgGradZ] = vectorGradient(filteredImg);
    imgGradSqNorm = sum(imgGradX.^2,4)+sum(imgGradY.^2,4)+sum(imgGradZ.^2,4);
    c = exp(-imgGradSqNorm/conductance);
    %c = c(:,:,:,ones(size(filteredImg,4),1));
    %P = vectorDivergence(c.*imgGradX, c.*imgGradY, c.*imgGradZ);
    P = weightedVectorDivergence(filteredImg, c)*1e2;
    
    step = 1;
    %g = B+priorWeight*P;
    g = P;
    %disp(['norm B: ', num2str(norm(B(:)))]);
    %disp(['norm P: ', num2str(norm(P(:)))]);
    %step = 1e-2;
    while step>1e-20 && logPosterior(filteredImg+step*g, noisyImg, sqStdDeviation, conductance, priorWeight)...
            <=logPosterior(filteredImg, noisyImg, sqStdDeviation, conductance, priorWeight)+a*step*dot(g(:),g(:))
        step = b*step;
    end
    disp(['step = ', num2str(step)]);
    filteredImg = filteredImg + step*g;
    
    logPost = logPosterior(filteredImg, noisyImg, stdDeviation, conductance, priorWeight);
    %disp(['posterior:', num2str(logPost)]);
    disp(['e:', num2str(logPost-lastLogPost)]);
    if abs(logPost-lastLogPost) <tolerance
        break;
    end;
    
    it = it+1;
    
end

%output
return;

end

function p = logPosterior(filteredImg, noisyImg, sqStdDeviation, conductance, priorWeight)

z = filteredImg.*noisyImg/sqStdDeviation;
p1 = log(noisyImg/sqStdDeviation)-(filteredImg.^2+noisyImg.^2)/2/sqStdDeviation+log(besseli(0,z,1))+abs(z);
p1 = sum(p1(:));

[imgGradX, imgGradY, imgGradZ] = vectorGradient(filteredImg);
imgGradSqNorm = sum(imgGradX.^2,4)+sum(imgGradY.^2,4)+sum(imgGradZ.^2,4);
E = exp(-imgGradSqNorm/conductance);
E = priorWeight*sum(E(:));
%p2 = 1/(1-exp(-priorWeight*numel(noisyImg))) * exp(-E);
p2 = E;

%p = p1+p2;
p=p2;

end

function div = weightedVectorDivergence (vol, weight)

ndim = length(size(vol))-1;

vol = permute(vol, [ndim+1, 1:ndim]);
weight = weight(:,:,:,ones(size(vol,1),1));
weight = permute(weight, [ndim+1, 1:ndim]);
div = zeros(size(vol));

if ndim == 1
  perm = [1 2];
else
  perm = [2:ndim 1]; % Cyclic permutation
end

for k=1:ndim
    [vecdim, n, p] = size(vol);
    d = zeros(size(vol), class(vol));
    
    if n>1
        d(:,1,:) = d(:,1,:) + weight(:,2,:).*(vol(:,2,:) - vol(:,1,:));        
        d(:,n,:) = d(:,n,:) + weight(:,n-1,:).*(vol(:,n-1,:) - vol(:,n,:));
    end
    
    if n>2
        d(:,2:n-1,:) = d(:,2:n-1,:) + weight(:,3:n,:).*(vol(:,3:n,:)-vol(:,2:n-1,:));
        d(:,2:n-1,:) = d(:,2:n-1,:) + weight(:,1:n-2,:).*(vol(:,1:n-2,:)-vol(:,2:n-1,:));
    end
    
    div = div+ipermute(d,[0 k:ndim 1:k-1]+1);
    
    vol = permute(vol, [1 perm+1]);
    weight = permute(weight, [1 perm+1]);
end

div = ipermute(div, [ndim+1, 1:ndim]);

end