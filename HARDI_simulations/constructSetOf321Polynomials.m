function C=constructSetOf321Polynomials(order)

UnitVectors;
Mprime=321;
g=g(1:Mprime,:);
for i=0:order
    for j=0:order-i
        pop(i+1,j+1,order-i-j+1)=population(i,j,order-i-j,order);
    end
end
for k=1:length(g)
    c=1;
    for i=0:order
		for j=0:order-i
			C(k,c)=pop(i+1,j+1,order-i-j+1)*(g(k,1)^i)*(g(k,2)^j)*(g(k,3)^(order-i-j));
			c=c+1;
        end
    end
end
