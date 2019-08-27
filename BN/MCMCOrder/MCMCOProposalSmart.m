function [neworder,deltalog,score_neworder,bayes_new]=MCMCOProposalSmart(SCORES,order,bayes_orig);
%%% given scores and an order, propose a new order and calculate the
%%% difference in the two scores. The difference is calculated only in
%%% interval between the nodes that are being fliped.

%% create new order... get the position where this occurs n1,n2
[neworder,n1,n2]=MCMCONewOrder(order);
%% order the nodes%%%
x=sort([n1,n2]);
x1=x(1);
x2=x(2);

[bayes_new]=MCMCOFlip2nodes(SCORES,neworder,order,n1,n2,bayes_orig);
score_order=MCMCOmatrix2score(bayes_orig);
score_neworder=MCMCOmatrix2score(bayes_new);
deltalog=score_neworder-score_order;

