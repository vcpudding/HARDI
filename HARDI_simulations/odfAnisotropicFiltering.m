function [filteredODF eBuf] = odfAnisotropicFiltering (noisyODF, conductance, tolerance, step)

filteredODF = noisyODF;
e = energyFunc(filteredODF, conductance);
eBuf =[];
it=1;

while it<=200
    laste = e;
    dODF = odfAnisotropicDivergence(filteredODF, conductance);
    filteredODF = expMap(filteredODF, step*dODF);
    e = energyFunc(filteredODF, conductance);   
   
    disp(['#', num2str(it), '  e:', num2str(e)]);
    eBuf = [eBuf e];
    
    if e-laste < tolerance
        break;
    end
    it = it+1;
end

end

function e = energyFunc (odf, conductance)

g = odfGradientSqNorm(odf);
e = exp(-g/conductance);
e = sum(e(:));

end

function l = logMap(odf1, odf2)

d = dot(odf1, odf2);
l = (odf2-d*odf1)/(1-d^2)*acos(d);

end

function e = expMap(odf, dODF)

d = sqrt(sum(odf.^2, 4));
d = d(:,:,:,ones(size(odf,4),1));
e = cos(d).*odf + sin(d).*dODF./d;

end

function g = odfGradientSqNorm (odf)

[m n p dim] = size(odf);
g = zeros(m,n,p);

for i=1:m
    for j=1:n
        for k=1:p
            if i==1
                gy = logMap(odf(i+1,j,k,:), odf(i,j,k,:));
            else
                if i==m
                    gy = logMap(odf(i,j,k,:), odf(i-1,j,k,:));
                else
                    gy = logMap(odf(i+1,j,k,:), odf(i-1,j,k,:))/2;
                end
            end
               
            if j==1
                gx = logMap(odf(i,j+1,k,:), odf(i,j,k,:));
            else
                if j==n
                    gx = logMap(odf(i,j,k,:), odf(i,j-1,k,:));
                else
                    gx = logMap(odf(i,j+1,k,:), odf(i,j-1,k,:))/2;
                end
            end
               
            if k==1
                gz = logMap(odf(i,j,k+1,:), odf(i,j,k,:));
            else
                if k==p
                    gz = logMap(odf(i,j,k,:), odf(i,j,k-1,:));
                else
                    gz = logMap(odf(i,j,k+1,:), odf(i,j,k-1,:))/2;
                end
            end
            
            g(i,j,k) = dot(gy,gy)+dot(gx,gx)+dot(gz,gz);
        end
    end
end

end

function div = odfAnisotropicDivergence (odf, conductance)

[m n p dim] = size(odf);
div = zeros(size(odf));

for i=1:m
    for j=1:n
        for k=1:p
            gN = logMap(odf(max(1,i-1),j,k,:), odf(i,j,k,:));
            gS = logMap(odf(min(m,i+1),j,k,:), odf(i,j,k,:));
            gW = logMap(odf(i,max(1,j-1),k,:), odf(i,j,k,:));
            gE = logMap(odf(i,min(n,j+1),k,:), odf(i,j,k,:));
            gF = logMap(odf(i,j,max(1,k-1),:), odf(i,j,k,:));
            gP = logMap(odf(i,j,min(p,k+1),:), odf(i,j,k,:));
            
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