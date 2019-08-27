function [NewBayes]=MCMCOFlip2nodes(SCORES,neworder,order,n1,n2,bayes)
%% calculate the scores only using the nodes that are between nodes that
%% have been fliped.
%% bayes is a matrix with all the likelihoods from the previous order.

%order the nodes
x=[n1 n2];
x=sort(x);
x1=x(1);
x2=x(2);
NewBayes=bayes;
% to the first node we need to remove all the contributions from families
% that have any parent between x1+1 till x2 and that have been valid in
% previous order


%%% Node in x1: check if the node has any parents between x1+1 and x2. If
%%% have and was valid, then remove their contribution.

[NewBayes]=MCMCOCheckx1(SCORES,NewBayes,neworder,order,x1,x2);
[NewBayes]=MCMCOCheckx2(SCORES,NewBayes,neworder,order,x1,x2);
[NewBayes]=MCMCOCheckx1_x2(SCORES,NewBayes,neworder,order,x1,x2);

