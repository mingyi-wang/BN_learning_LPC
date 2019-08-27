clear

%% Load the pre-computed SCORES structure.
load SCORES
rand('state',35);

%% Define the number of MCMC steps and sampling interval
NroMcmcSteps=10000;
SampleInterval=10;

%% sample orders
Results=MCMCORunSmart(SCORES,NroMcmcSteps,SampleInterval);

%% In Results you will find:
%% sampledSteps  
%% acceptRatio
%% logOrderScore
%% sampled_order

%% you can save these results if you want to:
%% save MyResults Results


%% Now we sample DAGs given sampled orders
%% You need to provide the SCORES structure, the sampled orders and a
%% number specifying how many DAGs you want to sample for each order in
%% this example we choose 20

[DAGs,logDAGscores]=OrderMCMC_SampleDAGsGivenOrder(SCORES,Results.sampled_order,20);

%% Now you have your trajectory of DAGs you can store it if you wish.

%save MyDAGsResults DAGs logDAGscores

