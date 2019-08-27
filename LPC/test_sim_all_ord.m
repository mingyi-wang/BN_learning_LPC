function test_sim_all_ord

s=[10 20 50 100 200 500 1000];
genes=1000;
ord=-1;
for j=1:10
  for i=1:7

    fname= sprintf('%s%d%s%d%s%d%s%d%s','..\data\data_',genes,'genes\',genes,'g_10rec_x_',j,'_',s(i),'measurements.txt');
    fname1= sprintf('%s%d%s%d%s%d%s','..\results\lpc\',genes,'g_10rec_x_',j,'_',s(i),'_lpc_ug_all_ord.mat');
    x=dlmread(fname,'\t');
    
    [undirected_G,sep,zMin,ci_num] = learn_struct_lpc_for_test(x, ord, 0.05);
%     [pdag,sep] = Orientation_lpc(x, undirected_G, sep, ord, 0, 0.05);
%     [ii,jj]=find(pdag==1);
%     for k=1:length(ii)
%        pdag(ii(k),jj(k))=2;
%     end
%     [ii,jj]=find(pdag==-1);
%     for k=1:length(ii)
%        pdag(ii(k),jj(k))=1;
%     end   

   save(fname1,'undirected_G','sep','zMin','ci_num');
    i
    j

end
end
end