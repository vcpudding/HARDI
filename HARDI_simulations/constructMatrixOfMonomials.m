function G=constructMatrixOfMonomials(g,order)

for k=1:length(g)
    c=1;
    for i=0:order
		for j=0:order-i
			G(k,c)=(g(k,1)^i)*(g(k,2)^j)*(g(k,3)^(order-i-j));
			c=c+1;
        end
    end
end
