clear 
% Get the Trajectory of DAGs
load ResultsExample

% Get the posterior probabilities for DGE and DGE
DAGmDGE=media(DAGs,1);
DAGmUGE=media_UGE(DAGs);

%% Load the true netowork, the network which the data was generated from
idxtm=load('true_mat_hptc.txt');

% Get the AUC values
DGEAurocValueNoDiag2(DAGmDGE,idxtm,1,1)
UGEAurocValueNoDiag2(DAGmUGE,idxtm,1,1)

% Get the ROC curves with the AUC values
DGERocCurveNoDiag2(DAGmDGE,idxtm,1,1)
UGERocCurveNoDiag2(DAGmUGE,idxtm,1,2)
