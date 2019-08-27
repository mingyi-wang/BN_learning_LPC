function [true_m]=DGETrueMatrix(size_m,idxtm)
%%% creating true matrix, matrix with was used to generate data.
%load a file with nodes connections

true_m=zeros(size_m);

%% fill ones where is indicated by idxtm
for j=1:size(idxtm,1)
    true_m(idxtm(j,1),idxtm(j,2))=1;
end
