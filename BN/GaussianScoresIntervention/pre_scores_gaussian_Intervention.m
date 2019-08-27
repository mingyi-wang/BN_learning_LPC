% ----------------
% Preliminaries
% ----------------
clear
rand('seed',28);
flag_prior=1; % restrictive prior

% ----------------
% Data
% ----------------


% Load here your data and the vector with interventions.
% The data should have nodes in rows and observations in columns.
% The vector with interventions should have as many entries as there are
% observations in the data set. Each entry is a number that indicates
% which node was intervened in each observation. If there is no node
% intervened the entry should be zero. For instance the vector [0 0 1 2 0
% 0]' will indicate that there was no interventions in observations 1,2, 5
% and 6 and that node 1 was intervened in observation 3 and node 2 was
% intervened in observation 3.

% The function below will calculate all the prior parameters with the less informative prior.
% If you want to change the parameters you need to edit the file Compute_Prior_Info.m


DATA = data_to_DATA(data,Intervention);


[T_0, T_m, v, alpha] = Compute_Prior_Info_Intervention(DATA,data);


% --------------------------------
% Precomputation of scores
% --------------------------------
Nbest=176;
minDifference=1400;
SCORES=SBmcmc_MakeScoreStructureGaussianIntervention(DATA,data,T_0,T_m,v,alpha,'Nbest',Nbest,'minDifference',minDifference,'flag_prior',flag_prior);

save YOUR_SCORES_FILE SCORES