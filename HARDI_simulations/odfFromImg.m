function odf = odfFromImg (img)

load GradientOrientations_64
UnitVectors
order = 4;
delta = 200;
G=constructMatrixOfMonomials(g,order);

dims = size(img);
dims = dims(1:3);
img = permute(img, [4 1 2 3]);
odf = zeros(length(G), prod(dims));

for i=1:prod(dims)            
    S = img(:,i);
    [ODF_coef, ODF_tensor] = Estimate_tensorODF(S, 1, GradientOrientations, order, delta);
    odf(:,i) = G*ODF_coef;
    odf(:,i) = odf(:,i)/sum(odf(:,i));
    disp(i);
end        

odf = reshape(odf, [length(G) dims]);
odf = ipermute(odf, [4 1 2 3]);

end