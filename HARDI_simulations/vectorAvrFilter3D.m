function filteredImg = vectorAvrFilter3D (img, maskSize)

[m n p vecdim] = size(img);
s = (maskSize-1)/2;

filteredImg = zeros(size(img));

for i=1:m
    for j=1:n
        for k=1:p
            dat = img(max(1, i-s):min(m,i+s),...
                max(1, j-s):min(n, j+s),...
                max(1, k-s):min(p, k+s), :);
            dat = reshape(permute(dat,[4,1,2,3]), [vecdim, numel(dat)/vecdim]);
            filteredImg(i,j,k, :) = mean(dat,2);
        end
    end
end

end