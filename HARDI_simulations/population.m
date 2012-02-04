function counter=population(i,j,k,order)

    size=3^order;
    counter=0;
    for z=0:size-1
        c=populationBasis(z,order,[0 0 0]);
        if (c(1)==i)&(c(2)==j)&(c(3)==k)
            counter=counter+1;
        end
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ret=populationBasis(i,order,c)
	if order==0  
        ret=c;
    else
	    j=mod(i,3);
	    c(j+1)=c(j+1)+1;
	    ret=populationBasis((i-j)/3,order-1,c);
    end