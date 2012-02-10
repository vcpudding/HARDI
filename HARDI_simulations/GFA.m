function y = GFA (D)
D = D/sum(D);
n=length(D);
r=norm(D)/sqrt(n);
y=std(D)/r;
end