function [resul]=media_UGE(dagtrajectory)

%%% This function calculates the mean posterior over the edges from a resulting MCMC
%%% trajectory of DAGs. Only the skeleton is taken in account all edge
%%% directions aren't used.

r=dagtrajectory;
acc=0;
for k=1:size(dagtrajectory,2)
    acc=acc+(r{k}+r{k}');
end
resul=acc/size(dagtrajectory,2);

