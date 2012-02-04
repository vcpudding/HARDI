function y = GFA (D)
n=length(D);
r=norm(D)/sqrt(n);
y=std(D)/r;
end