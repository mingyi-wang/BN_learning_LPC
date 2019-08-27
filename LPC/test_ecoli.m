function test_ecoli
% 
% Run the LPC-algorithm over the E.coli microarray dataset, which is downloaded from http://gardnerlab.bu.edu/data/PLoS_2007/
%
load ('..\data\ecoli\compendium_E_coli_v3_Build_1.mat');
load('../data/ecoli/tfs/tfs.mat');
data=compendium.rma;
data=data(compendium.gidx,:);
p=size(data,1);
G=zeros(p,p);
G(tfs.gidx_4345,:)=1;
G(:,tfs.gidx_4345)=1;


tmp2=cputime;
[undirected_G,sep,zMin] = learn_struct_lpc(data', 3, 0.001,G);
tmp2=cputime-tmp2;
save('..\data\ecoli\pdag_ecoli_udag_ord3','undirected_G','sep','zMin','tmp2');

[pdag,sep] = Orientation_lpc(data, undirected_G, sep, 3, 1, 0.001);
tmp2=cputime-tmp2;
save('..\data\ecoli\pdag_ecoli_dag_ord3','pdag','sep','tmp2');
end 