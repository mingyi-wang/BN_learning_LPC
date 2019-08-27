
s=[10 20 50 100 200 500 1000];

num_a = zeros(10,10);
num_edges=zeros(10,10);
num_true = zeros(10,10);
num_mo = zeros(10,10);
num_ea = zeros(10,10);
num_ma = zeros(10,10);
num_wo = zeros(10,10);
num_shd = zeros(10,10);
precision=zeros(10,10);
recall=zeros(10,10);
acc=zeros(10,10);
mean_num_mo=[];
   mean_num_ea=[];
   mean_num_ma=[];
   mean_num_wo=[];
   mean_num_shd=[];
   
   std_num_mo=[];
   std_num_ea=[];
   std_num_ma=[];
   std_num_wo=[];
   std_num_shd=[];

for i=1:7
    for j=1:10
    %fname2=sprintf('%s%d%s%d%s','H:\buf\NetRec\100g_10rec_scale-free_A_',j,'.mat');
    fname2=sprintf('%s%d%s','..\data\data_100genes\100g_10rec_A_',j,'.txt');
    A=dlmread(fname2,'\t');
    [ii,jj]=find(A~=0);
    for kkk=1:length(ii)
        if A(ii(kkk),jj(kkk))==-1
          A(ii(kkk),jj(kkk))=1;          
        end
    end

   fname1= sprintf('%s%d%s%d%s','..\results\lpc\100g_10rec_x_',j,'_',s(i),'_lpc_all_ord_pdags.mat');
   load(fname1);
   for k=1:10
       if length(pdag{k})==0
           pdag{k}=pdag{k-1};
       end
   [num_a(j,k),num_edges(j,k),num_true(j,k),num_mo(j,k),num_ea(j,k),num_ma(j,k),num_wo(j,k),num_shd(j,k),recall(j,k),precision(j,k),acc(j,k)]=compare_g(pdag{k},A);
   end
   end
   mean_num_mo=[mean_num_mo; mean(num_mo,1)];
   mean_num_ea=[mean_num_ea; mean(num_ea,1)];
   mean_num_ma=[mean_num_ma; mean(num_ma,1)];
   mean_num_wo=[mean_num_wo; mean(num_wo,1)];
   mean_num_shd=[mean_num_shd; mean(num_shd,1)];
   
   std_num_mo=[std_num_mo;std(num_mo,0,1)];
   std_num_ea=[std_num_ea;std(num_ea,0,1)];
   std_num_ma=[std_num_ma;std(num_ma,0,1)];
   std_num_wo=[std_num_wo;std(num_wo,0,1)];
   std_num_shd=[std_num_shd;std(num_shd,0,1)];
 
   i
   j
end

figure
for i=1:7
  subplot(2,2,1);
  hold on;
  switch(i)
      case 1
         errorbar([1:10], mean_num_shd(i,:), std_num_shd(i,:), 'color',[1 0 1]);
      case 2
         errorbar([1:10], mean_num_shd(i,:), std_num_shd(i,:), 'color',[0 1 1]);
      case 3
         errorbar([1:10], mean_num_shd(i,:), std_num_shd(i,:), 'color',[1 0 0]);
      case 4
         errorbar([1:10], mean_num_shd(i,:), std_num_shd(i,:), 'color',[0 1 0]);
      case 5
         errorbar([1:10], mean_num_shd(i,:), std_num_shd(i,:), 'color',[0 0 1]);
      case 6
         errorbar([1:10], mean_num_shd(i,:), std_num_shd(i,:), 'color',[0 0.5 0.3]);
      case 7
         errorbar([1:10], mean_num_shd(i,:), std_num_shd(i,:), 'color',[0 0 0]);
  end
  subplot(2,2,2);
  hold on;
  switch(i)
      case 1
         errorbar([1:10], mean_num_mo(i,:), std_num_mo(i,:), 'color',[1 0 1]);
      case 2
         errorbar([1:10], mean_num_mo(i,:), std_num_mo(i,:), 'color',[0 1 1]);
      case 3
         errorbar([1:10], mean_num_mo(i,:), std_num_mo(i,:), 'color',[1 0 0]);
      case 4
         errorbar([1:10], mean_num_mo(i,:), std_num_mo(i,:), 'color',[0 1 0]);
      case 5
         errorbar([1:10], mean_num_mo(i,:), std_num_mo(i,:), 'color',[0 0 1]);
      case 6
         errorbar([1:10], mean_num_mo(i,:), std_num_mo(i,:), 'color',[0 0.5 0.3]);
      case 7
         errorbar([1:10], mean_num_mo(i,:), std_num_mo(i,:), 'color',[0 0 0]);
  end
  subplot(2,2,3);
  hold on;
  switch(i)
      case 1
         errorbar([1:10], mean_num_ea(i,:), std_num_ea(i,:), 'color',[1 0 1]);
      case 2
         errorbar([1:10], mean_num_ea(i,:), std_num_ea(i,:), 'color',[0 1 1]);
      case 3
         errorbar([1:10], mean_num_ea(i,:), std_num_ea(i,:), 'color',[1 0 0]);
      case 4
         errorbar([1:10], mean_num_ea(i,:), std_num_ea(i,:), 'color',[0 1 0]);
      case 5
         errorbar([1:10], mean_num_ea(i,:), std_num_ea(i,:), 'color',[0 0 1]);
      case 6
         errorbar([1:10], mean_num_ea(i,:), std_num_ea(i,:), 'color',[0 0.5 0.3]);
      case 7
         errorbar([1:10], mean_num_ea(i,:), std_num_ea(i,:), 'color',[0 0 0]);
  end
  subplot(2,2,4);
  hold on;
  switch(i)
      case 1
         errorbar([1:10], mean_num_ma(i,:), std_num_ma(i,:), 'color',[1 0 1]);
      case 2
         errorbar([1:10], mean_num_ma(i,:), std_num_ma(i,:), 'color',[0 1 1]);
      case 3
         errorbar([1:10], mean_num_ma(i,:), std_num_ma(i,:), 'color',[1 0 0]);
      case 4
         errorbar([1:10], mean_num_ma(i,:), std_num_ma(i,:), 'color',[0 1 0]);
      case 5
         errorbar([1:10], mean_num_ma(i,:), std_num_ma(i,:), 'color',[0 0 1]);
      case 6
         errorbar([1:10], mean_num_ma(i,:), std_num_ma(i,:), 'color',[0 0.5 0.3]);
      case 7
         errorbar([1:10], mean_num_ma(i,:), std_num_ma(i,:), 'color',[0 0 0]);
  end
end
  
