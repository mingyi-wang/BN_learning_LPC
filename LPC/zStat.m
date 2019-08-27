function [res,r]=zStat(x,y,S,C,n)
r=pcorOrder(x,y,S,C);
res=sqrt(n-length(S)-3).*(0.5.*log((1+r)./(1-r)));
if isinf(res)
    res=0;
end
end