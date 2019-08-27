
% ----------------
% Data for Simulated datasets from Italy programs
% ----------------

% Load here your data.

type='scale-free';
%type='random';
%type='scale-free_timecourse';
s=[10 20 50 100 200 500 1000];
% 
%   data1=dlmread('f:\projects\e_tfs_nogene.csv',',');
%   data2=dlmread('f:\projects\e_transporter_nogene.csv',',');
%   
%   data=[data1;data2];
%   data=data';
%fid=fopen('h:\buf\NetRec\scale_free_results.txt','w+');
% fid=fopen('h:\buf\Software_Germany\BN_results.txt','w+');
for j=2:2
    %fname2=sprintf('%s%d%s%d%s','H:\buf\NetRec\100g_10rec_scale-free_A_',j,'.mat');
%     fname2=sprintf('%s%s%s%d%s','H:\buf\NetRec\data\100g_10rec_',type,'_A_',j,'.txt');
%     A=dlmread(fname2,'\t');
    
for i=1:1
  %  if (j==3 & i==1) |(j==3 & i==3) | ((j==4) & (i==1|i==2|i==3) )
        
            
%    if i<10 
%       fname= strcat('H:\projects\constraint_bn\data\AGN-Century\CenturyRND-00',num2str(i),'.nmss.dat');
%       fname1= strcat('H:\projects\constraint_bn\data\AGN-Century\CenturyRND-00',num2str(i));
%       fname='';
%       fname1='';
%       fname= sprintf('%s%d%s%d%s','H:\buf\NetRec\100g_10rec_scale-free_x_',j,'_',s(i),'measurements.mat');
%       fname1= sprintf('%s%d%s%d%s','H:\buf\NetRec\100g_10rec_scale-free_x_',j,'_',s(i),'_lowPC.mat');
      fname= sprintf('%s%d%s%d%s','H:\projects\LPC\data\data_100genes\100g_10rec_x_',j,'_',s(i),'measurements.txt');
      fname1= sprintf('%s%d%s%d%s','H:\projects\LPC\results\bn\100g_10rec_x_',j,'_',s(i),'_bn.mat');

%     else
%       fname= strcat('H:\projects\constraint_bn\data\AGN-Century\CenturyRND-0',num2str(i),'.nmss.dat');
%       fname1= strcat('H:\projects\constraint_bn\data\AGN-Century\CenturyRND-0',num2str(i));  
%     end
    x=dlmread(fname,'\t');
    x=x';
%It should have nodes in rows and observations in columns
%In this case we are going to use the first data set that was generated
%from a multivarite gaussian with observations only.


% The function below will calculate all the prior parameters with the less informative prior.
% if you want to change the parameters you need to edit the file Compute_Prior_Info.m
cd('H:\projects\LPC\bn\GaussianScores');
%x=x(1:10,:);
[T_0, T_m, v, alpha] = Compute_Prior_Info(x);


% --------------------------------
% Precomputation of scores
% --------------------------------

% The scores are calculated for possible families for each node. For more information type help SBmcmc_MakeScoreStructureGaussian

Nbest=176;
minDifference=1000;
flag_prior=1;
%SCORES=SBmcmc_MakeScoreStructureGaussian(x,T_0,T_m,v,alpha,'Nbest',Nbest,'minDifference',minDifference,'flag_prior',flag_prior);
SCORES=SBmcmc_MakeScoreStructureGaussian(x,T_0,T_m,v,alpha,'flag_prior',flag_prior);
cd('H:\projects\LPC\bn\MCMCOrder');
NroMcmcSteps=10000;
SampleInterval=10;

Results=MCMCORunSmart(SCORES,NroMcmcSteps,SampleInterval);


[DAGs,logDAGscores]=OrderMCMC_SampleDAGsGivenOrder(SCORES,Results.sampled_order,20);

%% Now you have your trajectory of DAGs you can store it if you wish.

save (fname1,'SCORES','DAGs','logDAGscores');
clear SCORES
clear DAGs
clear logDAGscores
        
    end
%end
end

