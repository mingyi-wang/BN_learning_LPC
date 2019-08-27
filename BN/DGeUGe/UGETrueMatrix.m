function [true_m]=UGETrueMatrix(size_m,idxtm)
%%% creating true matrix, matrix with was used to generate data.
%load a file with nodes connections

true_m=zeros(size_m);% create a matrix with diagonal 1

%% fill ones where is indicated by idxtm and also in the symmetric
%% positions
for j=1:size(idxtm,1)
    true_m(idxtm(j,1),idxtm(j,2))=1;
    true_m(idxtm(j,2),idxtm(j,1))=1;
end
