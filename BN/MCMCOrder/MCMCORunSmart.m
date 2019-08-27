function [Results]=MCMCORunSmart(SCORES,NroSteps,SampleInterval)


nro_nodes=length(SCORES);
order=randperm(nro_nodes);
[bayes_orig]=MCMCOTotLog(SCORES,order);
nro_proposal=0;
nro_accepted=0;
nro_sample=0;
nro_sampled_order=0;


for step=1:NroSteps
    %step
    nro_proposal=nro_proposal+1;
    [neworder,deltalog,score_neworder,bayes_new]=MCMCOProposalSmart(SCORES,order,bayes_orig);
    if deltalog>=0
        flag_accepted=1;
    elseif exp(deltalog)>rand
        flag_accepted=1;
    else
        flag_accepted=0;
    end
    
    if flag_accepted==1
        nro_accepted=nro_accepted+1;
        order=neworder;
        bayes_orig=bayes_new;
    end
    
    %% keeping results
    if step==1 | mod(step,SampleInterval)==0
     nro_sample=nro_sample+1;  
     sampledSteps(nro_sample)=step;
     acceptRatio(nro_sample)=nro_accepted/nro_proposal;
     nro_proposal=0;
     nro_accepted=0;
     logOrderScore(nro_sample)=score_neworder;
     if step>NroSteps/2
        nro_sampled_order=nro_sampled_order+1;
        sampled_order(nro_sampled_order,:)=neworder;
     end
    end
end%mcmc

Results.sampledSteps=sampledSteps;
Results.acceptRatio= acceptRatio;
Results.logOrderScore=logOrderScore;
Results.sampled_order=sampled_order;
        
    