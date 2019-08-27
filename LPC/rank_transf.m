% rank transformation: returns vector of integers corresponding to their
% ordering in the original input vector
function ib = rank_transf(x)
[a ia] = sort(x);
[b ib] = sort(ia);
