% Calcuate the mean CI test numbers and the SHD values at order=0-9 
%
ss=[10 20 50 100 200 500 1000];  % Sample sizes

mean_cis=[];
std_cis=[];
mean_cis2=[];
std_cis2=[];
mean_aucs=[];
std_aucs=[];

max_ord=zeros(7,10);
hold all;
for i=1:7 % sample size
  ci=zeros(20,10);
  ci_2=zeros(3,10);
  for j=1:10  % repeat tests
    fname=['../data/data_100genes/100g_10rec_A_', num2str(j),'.txt'];
    A=dlmread(fname,'\t');
    A_undir=A|A';
    filename= ['../results/lpc/100g_10rec_x_' num2str(j), '_',num2str(ss(i)),'_lpc_ug_all_ord.mat'];
    load (filename);
    for k=1:10  % order
      if ci_num(k)==0
        break; 
      end
    end
    max_ord(i,j)=k-1;
    for k=1:max_ord(i,j)
      k
      if k>1 
        ci(k,j)=ci(k-1,j)+ci_num(k);
      else
        ci(k,j)=ci_num(k);         
      end
    end
    for k=1:max_ord(i,j)
      k
      if k>3
        ci_2(2,j)=ci_2(2,j)+ci_num(k);
      else
        ci_2(1,j)=ci_2(1,j)+ci_num(k); 
      end
      ci_2(3,j)=ci_2(3,j)+ci_num(k);
    end
  end
  mean_ci=mean(ci,2);
  std_ci=std(ci,0,2);
  mean_ci2=mean(ci_2,2);
  std_ci2=std(ci_2,0,2);
  
  mean_cis=[mean_cis,mean_ci];
  std_cis=[std_cis,std_ci];
  
  mean_cis2=[mean_cis2,mean_ci2];
  std_cis2=[std_cis2,std_ci2];
  hold on; 
   
end
