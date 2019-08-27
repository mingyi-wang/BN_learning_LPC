function [resul]=media(dagtrajectory,pdag)
%%% function [resul]=media(dagtrajectory,pdag)
%%% 
%%% This function calculates the mean posterior over the edges from a resulting MCMC
%%% trajectory of DAGs. There is an option to convert to PDAG or not...
%%% dagtrajectory.
%%% A collection of DAGs sampled during an MCMC simulation.
%%%
%%% pdag
%%% if pdag is 1, each DAG from the trajectory is converted to a pdag
%%% before proceeding to the mean calculation. Any other number will cause
%%% the algorithm to not convert to PDAGs

if pdag==1
    sprintf('I am converting DAGs for PDAG...')
    why
    for m=1:size(dagtrajectory,2)
        r{m}=ConversionDAG2PDAG(dagtrajectory{m});
        if mod(m,100)==0
            m
        end
    end
else
    sprintf('I am not converting DAGs for PDAG...')
    why
    r=dagtrajectory;
end

acc=0;
for k=1:size(dagtrajectory,2)
    acc=acc+r{k};
end
resul=acc/size(dagtrajectory,2);

