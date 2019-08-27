clear 
% Get the Trajectory of DAGs
s=[10 20 50 100 200 500 1000];
i=1;  %Graph
type='scale-free';
% fname1= sprintf('%s%s%s%d%s%d','H:\buf\NetRec\data\bn\100g_10rec_',type,'_x_',j,'_',i,'_bn.mat');
% load(fname1);
% 
% % Get the posterior probabilities for DGE and DGE
% DAGmDGE=media(DAGs,1);
% DAGmUGE=media_UGE(DAGs);
% Get the AUC values
% DGEAurocValueNoDiag2(DAGmDGE,idxtm,1,1)
% UGEAurocValueNoDiag2(DAGmUGE,idxtm,1,1)
% 
% % Get the ROC curves with the AUC values
% DGERocCurveNoDiag2(DAGmDGE,idxtm,1,1)
% UGERocCurveNoDiag2(DAGmUGE,idxtm,1,2)
auc_rn=zeros(7,10);
auc_lowpc=zeros(7,10);
%% Load the true netowork, the network which the data was generated from
for i=1:10
fname=sprintf('%s%s%s%s%s','H:\buf\NetRec\data\100g_10rec_',type,'_A_',num2str(i),'.edg');
for j=1:7
idxtm=dlmread(fname,'\t');
output_ggm_directory = 'H:\buf\NetRec\data\ggm\';
output_lowpc_directory = 'H:\buf\NetRec\data\lowpc\';
output_rn_directory = 'H:\buf\NetRec\data\rn\';
filenames_start = [num2str(100) 'g_' num2str(10) 'rec_' type];
% 
% w_ggm=dlmread([output_ggm_directory filenames_start '_x_' num2str(i) '_' num2str(j) '_ggm.txt'],'\t');
% w_ggm=abs(w_ggm);
% DGEAurocValueNoDiag2(w_ggm,idxtm,1,1)
% UGEAurocValueNoDiag2(w_ggm,idxtm,1,1)
% 
% % Get the ROC curves with the AUC values
% DGERocCurveNoDiag2(w_ggm,idxtm,1,1)
% UGERocCurveNoDiag2(w_ggm,idxtm,1,2)

load([output_lowpc_directory filenames_start '_x_' num2str(i) '_' num2str(s(j)) '_lowPC_oldcc_ord2.mat']);
w_lowpc_oldcc=zMin;
%DGEAurocValueNoDiag2(w_lowpc_oldcc,idxtm,1,1)
auc_lowpc(j,i)=UGEAurocValueNoDiag2(w_lowpc_oldcc,idxtm,1,1)

% Get the ROC curves with the AUC values
%DGERocCurveNoDiag2(w_lowpc_oldcc,idxtm,1,1)
%UGERocCurveNoDiag2(w_lowpc_oldcc,idxtm,1,2)

load([output_rn_directory filenames_start '_x_' num2str(i) '_' num2str(s(j)) '_rn.mat']);
w_rn=RelPe;
for k=1:100
    w_rn(k,k)=0;
end
%DGEAurocValueNoDiag2(w_rn,idxtm,1,1)
auc_rn(j,i)=UGEAurocValueNoDiag2(w_rn,idxtm,1,1)

% Get the ROC curves with the AUC values
%DGERocCurveNoDiag2(w_rn,idxtm,1,1)
%UGERocCurveNoDiag2(w_rn,idxtm,1,2)
end
end
figure
plot(s,mean(auc_lowpc,2),'-+');
hold on;
plot(s,mean(auc_rn,2),'g');