function [r]=pcorOrder(i,j,k,C)
cutoff=0.9999999;
ord=length(k);
if (ord==0)    
    r=C(i,j);
elseif (ord==1)
    r=(C(i,j)-C(i,k).*C(j,k))./sqrt((1-C(j,k)^2).*(1-C(i,k)^2));
    
else
        s=k(ord);
        k=k(setdiff(1:ord,ord));
        w1=pcorOrder(i,j,k,C);
        w2=pcorOrder(i,s,k,C);
        w3=pcorOrder(j,s,k,C);
        r=(w1-w2.*w3)./sqrt((1-w2^2).*(1-w3^2));
end
if isinf(r)
    r=0;
end
r=min(cutoff,max(-cutoff,r));
end
 
    