function S = removeRicianNoiseSingleVoxel (S_noisy, sig, tol)

%S_noisy: MXN matrix, where M is the number of repititions and N is the
%number of gradients

[nGradients, nReps] = size(S_noisy);
sqSig = sig^2;

S = mean(S_noisy, 2);
a = 0.1;
b = 0.5;
it = 1;
%resBuf = [];
l = 0;
while 1
    S1 = S*ones(1, nReps);
    x = S1.*S_noisy/sqSig;
    B = (-S1/sqSig + (besseli(1,x,1)./besseli(0,x,1)).*S_noisy/sqSig)*ones(nReps, 1)/nReps;
    step = 1;
    while step>1e-20 && likelihoodFunc(S+step*B, S_noisy, sig) <= likelihoodFunc(S, S_noisy, sig)+a*step*dot(B,B)
        step = b*step;  
    end
    S = S+step*B;
    lastL = l;
    l = likelihoodFunc(S, S_noisy, sig);
    res = abs(lastL-l);
    if it>1 && res<tol
        break;
    end;
    %resBuf(it) = res;
    it = it+1;  
end

%figure, plot(log(resBuf));

end

function l = likelihoodFunc (S, S_noisy, sig)
[nGradients, nReps] = size(S_noisy);
sqSig = sig^2;
S1 = S*ones(1, nReps);
l = (log(S_noisy/sqSig)-(S_noisy.^2+S1.^2)/(2*sqSig)+log(besseli(0, S_noisy.*S1/sqSig, 1))+S_noisy.*S1/sqSig)*ones(nReps,1);
%l = (log(S_noisy/sqSig)-(S_noisy.^2+S1.^2)/(2*sqSig)+log(besseli(0, S_noisy.*S1/sqSig)))*ones(nReps,1);
l = mean(l);
end