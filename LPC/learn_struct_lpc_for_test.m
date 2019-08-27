function [undirected_G,sep_all,zMin_all,ci_num] = learn_struct_lpc_for_test(Data, k, alpha)
% LEARN_STRUCT_LPC Learn a undirected graph pattern using the low-order PC algorithm
% This algorithm can record the undirected graphs, sep sets, zMin, CI test numbers under different k, which is varied from 0 to highest order
% This algorithm is used for tests the "effect of the order on the LPC-algorithm
% 
% Data is matrix with n*p, n rows means observations and p means the nodes 
% k is an maximal order of PC algorithm (default: 2), k=-1 means that no limit on order of conditional independence tests
% alpha is the significant level for partial correlation, the default alpha value is 0.05
%
% The low order PC algorithm does structure learning assuming all variables are observed.
% See Spirtes, Glymour and Scheines, "Causation, Prediction and Search", 1993, p117.
% 
%


if nargin < 3, cutoff=norminv(1-0.05/2,0,1);
else
cutoff=norminv(1-alpha/2,0,1);
end  
if nargin <2, k=2; end % The highest order for CI tests

sep_all=cell(10,1);
zMin_all=cell(10,1);
undirected_G=cell(10,1);

p=size(Data,2);
n=size(Data,1);  
sep = cell(p,p);
ord = 0;

done = 0;
G = ones(p,p);
G=setdiag(G,0);

if k==-1
  ci_num=zeros(p-1,1);
else
  ci_num=zeros(k+1,1);
end
% This is just for ecoli data

% tf_ecoli=dlmread('./data/ecoli/ecoli_genes.txt','\t');
% gene_sn=find(tf_ecoli(:,2)==0)
% tf_sn=find(tf_ecoli(:,2)==1)
% for i=1:length(G)
%     if ismember(i, gene_sn)
%       G(i,gene_sn)=0;           
%     end
% end
    
%done=0;
zMin=inf(p,p);  %Fisher's z-transform of the partial correlation values for each node pair
zMin=setdiag(zMin,0);
% rMin=inf(p,p);
% rMin=setdiag(rMin,0);

C=corrcoef(Data);
while ~done
  done = 1;
  [X,Y] = find(G==1); 
  for i=1:length(X)
    x = X(i); y = Y(i);

    if x<y      
      nbrsx=find(G(x,:)>0);
      nbrs=nbrsx;
      nbrs = mysetdiff(nbrs, x);  % bug fix by Raanan Yehezkel <raanany@ee.bgu.ac.il> 6/27/04
      nbrs=mysetdiff(nbrs,y);
      if length(nbrs) >= ord & G(x,y) ~= 0
 
        SS = subsets1(nbrs, ord);% all subsets of size ord
        for si=1:length(SS)
	      S = SS{si};
          [z,r]=zStat(x,y,S,C,n);
          ci_num(ord+1)=ci_num(ord+1)+1;  % to record the number of CI tests performed
          %[rho,p]=partialcorr(Data(:,x),Data(:,y),Data(:,S));
          if abs(z)<zMin(x,y)
            zMin(x,y)=abs(z);
            zMin(y,x)=abs(z);
          end
          

          if abs(z)<=cutoff % | r<0   % This is change for medicago data analysis
             G(x,y) = 0;
	         G(y,x) = 0;
             sep{x,y} = S;
	       %  sep{y,x} = S;
            % done1=1;
             fprintf('%d indep of %d given ', x, y); fprintf('%d ', S);
             fprintf('\n');
             break;
          end         
        end
        if G(x,y)~=0 & length(nbrs) > ord  % for this pair x and y with specified order, they are still connected
          done=0;    % So higher ord CI tests needed
        end
      end       
    end    
  end
  ord = ord + 1;
  undirected_G{ord,1}=G;
  sep_all{ord,1}=sep;
  zMin_all{ord,1}=zMin;
  if ord>k & k~=-1
    done=1;
  end
  fprintf('order=%d\n',ord);
end

%undirected_G=G;


end


      