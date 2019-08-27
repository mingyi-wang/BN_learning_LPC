% ----------------
% Preliminaries
% ----------------
clear
rand('seed',17);
flag_prior=1; % restrictive prior

% ----------------
% Data
% ----------------

% Load here your data. It should have nodes in rows and observations in columns

% The function below will calculate all the prior parameters with the less informative prior.
% if you want to change the parameters you need to edit the file Compute_Prior_Info.m

[T_0, T_m, v, alpha] = Compute_Prior_Info(data);


% --------------------------------
% Precomputation of scores
% --------------------------------

% The scores are calculated for possible families for each node. For more information type help SBmcmc_MakeScoreStructureGaussian

Nbest=176;
minDifference=1000;
SCORES=SBmcmc_MakeScoreStructureGaussian(data,T_0,T_m,v,alpha,'Nbest',Nbest,'minDifference',minDifference,'flag_prior',flag_prior);


save YOUR_SCORES_FILE SCORES