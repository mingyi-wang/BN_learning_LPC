clear 
% Get the Trajectory of DAGs
AUC_DGEs=zeros(1,10);
AUC_UGEs=zeros(1,10);
recalls=zeros(1,10);
precisions=zeros(1,10);
accs=zeros(1,10);
SHDs=zeros(1,10);
fname2=sprintf('%s%d%s','H:\buf\NetRec\data\bn_data\bn_AUC.mat');

for j=1:10


fname1=sprintf('%s%d%s','H:\buf\NetRec\data\bn_data\bn_test_',j,'.mat');
load(fname1);
% Get the posterior probabilities for DGE and DGE
DAGmDGE=media(DAGs,1);
DAGmUGE=media_UGE(DAGs);
%% Load the true netowork, the network which the data was generated from

% idxtm=load('true_mat_hptc.txt');
% % Get the AUC values
% AUC_DGEs(j)=DGEAurocValueNoDiag2(DAGmDGE,idxtm,1,1)
% AUC_UGEs(j)=UGEAurocValueNoDiag2(DAGmUGE,idxtm,1,1)
% 
% % Get the ROC curves with the AUC values
% DGERocCurveNoDiag2(DAGmDGE,idxtm,1,1)
% UGERocCurveNoDiag2(DAGmUGE,idxtm,1,2)
% 
 
% dag=zeros(11,11);
%  for ii=1:length(idxtm)
%       dag(idxtm(ii,1),idxtm(ii,2))=1;
%       dag(idxtm(ii,2),idxtm(ii,1))=1;
%  end
 DAGmUGE(find(DAGmUGE>0.5))=1;
 
 DAGmUGE(find(DAGmUGE<=0.5))=0;
 DAGmDGE(find(DAGmDGE>0.5))=1;
 DAGmDGE(find(DAGmDGE<=0.5))=0;
 dag=dlmread('H:\buf\NetRec\data\sim_dat\A.txt','\t');
 dag(find(dag>0))=1;
 [ii,jj]=find(dag>0);
  idxtm=zeros(length(ii),2);
 for iiii=1:length(ii)
     idxtm(iiii,1)=ii(iiii);
     idxtm(iiii,2)=jj(iiii);
 end
 AUC_DGEs(j)=DGEAurocValueNoDiag2(DAGmDGE,idxtm,1,1)
 AUC_UGEs(j)=UGEAurocValueNoDiag2(DAGmUGE,idxtm,1,1)
 DGERocCurveNoDiag2(DAGmDGE,idxtm,1,1);
 UGERocCurveNoDiag2(DAGmUGE,idxtm,1,2);
cd ('H:\projects\constraint_bn');
[recall,precision,acc]=compare_uge(DAGmUGE,dag|dag');
recalls(j)=recall;
precisions(j)=precision;
accs(j)=acc;
[num_a,num_a1,num_true,num_mo,num_ea,num_ma,num_wo,SHD,recall,precision,acc]=compare_g(DAGmDGE,dag)
cd ('H:\buf\Software_Germany\DGeUGe\');
SHDs(j)=SHD;
end
save(fname2,'AUC_DGEs','AUC_UGEs','recalls','precisions','accs','SHDs');