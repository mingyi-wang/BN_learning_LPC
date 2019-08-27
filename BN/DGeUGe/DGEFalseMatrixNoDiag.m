function [false_m]=DGEFalseMatrixNoDiag(size_m,idxtm)
%%% creating false matrix, matrix with ones where there is no connections.
%load a file with nodes connections that represent the true connections

false_m=ones(size_m);% create a matrix with ones.

%% fill zeros where is indicated by idxtm
for j=1:size(idxtm,1)
    false_m(idxtm(j,1),idxtm(j,2))=0;
end
false_m=false_m-eye(size_m);