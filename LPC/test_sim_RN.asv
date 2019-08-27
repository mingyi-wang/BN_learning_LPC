function test_sim_RN


%type='random';
%type='small-world';
%type='scale-free_timecourse';
%type='random_timecourse';
s=[10 20 50 100 200 500 1000];
%s=[10 30 60 150 900];
% 
%   data1=dlmread('f:\projects\e_tfs_nogene.csv',',');
%   data2=dlmread('f:\projects\e_transporter_nogene.csv',',');
%   
%   data=[data1;data2];
%   data=data';
%fid=fopen('h:\buf\NetRec\scale_free_results.txt','w+');
% fname3=sprintf('%s%s%s','h:\buf\NetRec\data\lowpc\',type,'_results.txt');
% fid=fopen(fname3,'w+');
for j=1:10
    %fname2=sprintf('%s%d%s%d%s','H:\buf\NetRec\100g_10rec_scale-free_A_',j,'.mat');
%     fname2=sprintf('%s%s%s%d%s','H:\buf\NetRec\data\100g_10rec_',type,'_A_',j,'.txt');
%     A=dlmread(fname2,'\t');
for i=1:7
%    if i<10 
%       fname= strcat('H:\projects\constraint_bn\data\AGN-Century\CenturyRND-00',num2str(i),'.nmss.dat');
%       fname1= strcat('H:\projects\constraint_bn\data\AGN-Century\CenturyRND-00',num2str(i));
%       fname='';
%       fname1='';
%       fname= sprintf('%s%d%s%d%s','H:\buf\NetRec\100g_10rec_scale-free_x_',j,'_',s(i),'measurements.mat');
%       fname1= sprintf('%s%d%s%d%s','H:\buf\NetRec\100g_10rec_scale-free_x_',j,'_',s(i),'_lowPC.mat');
%      fname= sprintf('%s%s%s%d%s%d','H:\buf\NetRec\data\data_degree_50\100g_10rec_',type,'_x_',j,'_',s(i),'measurements.txt');
      fname= sprintf('%s%d%s%d%s','..\data\data_100genes\100g_10rec_x_',j,'_',s(i),'measurements.txt');
      fname1= sprintf('%s%d%s%d%s','..\results\rn\100g_10rec_x_',j,'_',s(i),'_RN.mat');

%     else
%       fname= strcat('H:\projects\constraint_bn\data\AGN-Century\CenturyRND-0',num2str(i),'.nmss.dat');
%       fname1= strcat('H:\projects\constraint_bn\data\AGN-Century\CenturyRND-0',num2str(i));  
%     end
    x=dlmread(fname,'\t');
    %data=dlmread(fname,'\t');
   RelPe=abs(corrcoef(x));
   
   save(fname1,'RelPe');
    i
    j
   
end
end
end