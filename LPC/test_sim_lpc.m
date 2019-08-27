function test_sim_lpc

s=[10 20 50 100 200 500 1000];
genes=100;
ord=2;
for j=1:1
  for i=1:7

    fname= ['..\data\data_',num2str(genes),'genes\',num2str(genes),'g_10rec_x_',num2str(j),'_',num2str(s(i)),'measurements.txt'];
    fname1= ['..\results\lpc\',num2str(genes),'g_10rec_x_',num2str(j),'_',num2str(s(i)),'_lpc_ord',num2str(ord),'.mat'];
    x=dlmread(fname,'\t');
    
    [undirected_G,sep,zMin,ci_num] = learn_struct_lpc(x, ord, 0.05);
    [pdag,sep] = Orientation_lpc(x, undirected_G, sep, ord,1, 0.05);
    
    save(fname1,'undirected_G','pdag','sep','zMin','ci_num');
    i
    j

end
end
end