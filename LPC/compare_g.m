function [num_a,num_a1,num_true,num_mo,num_ea,num_ma,num_wo,SHD,recall,precision,acc]=compare_g(G,dag)
% G is the DAG which want to be compared, G(x,y)=1 means x->y; G(x,y)=2 and
% G(y,x)=2 means x-y
% dag is the CPDAG
% 
% Structure of Hamming distance calculation
%
num_a=size(find(dag==1),1);


[II JJ]=find(G==1);
num_true=0;   %True prediction
for i=1:length(II)
  if (dag(II(i),JJ(i))==1) | (dag(JJ(i),II(i))==1)
    num_true=num_true+1;
  end
end
num_a1=length(II);


num_mo=0;     %Missing Orientation
[II JJ]=find(G==2);
for i=1:length(II)
  if II(i)<JJ(i)
    if ((dag(II(i),JJ(i))==1) & (dag(JJ(i),II(i))==0)) | ((dag(JJ(i),II(i))==1) & (dag(II(i),JJ(i))==0))
       
         num_mo=num_mo+1;
       
    end
  end
end
num_a1=num_a1+length(II)/2;

%Extra Arc
[II JJ]=find(G==1|G==2);
num_ea=0;
for i=1:length(II)
  if G(II(i),JJ(i))==2  % having no direction
    if II(i)<JJ(i)
       if dag(II(i),JJ(i))==0 & dag(JJ(i),II(i))==0
          num_ea=num_ea+1;
       end
    end
  else
     if dag(II(i),JJ(i))==0 & dag(JJ(i),II(i))==0
        num_ea=num_ea+1;    
     end
  end
end

%Missing Arcs
num_ma=0;
[II JJ]=find(dag==1);
for i=1:length(II)
    if dag(JJ(i),II(i))==1
       if II(i)<JJ(i)
         if G(II(i),JJ(i))==0 & G(JJ(i),II(i))==0
           num_ma=num_ma+1;
         end
       end 
    else
       if G(II(i),JJ(i))==0 & G(JJ(i),II(i))==0
         num_ma=num_ma+1;
       end        
    end
end


%Wrongly oriented
num_wo=0;
[II JJ]=find(dag==1);
for i=1:length(II)
   if dag(JJ(i),II(i))~=0
    if G(JJ(i),II(i))==1
         num_wo=num_wo+1;
    end
   end
end
SHD=num_wo+num_ma+num_ea+num_mo;
%fprintf('num:%d, true:%d, M.O.:%d, E.A.:%d, M.A.:%d, W.O.:%d\n',num_a,num_true,num_mo,num_ea,num_ma,num_wo);


%recall for undirected graph
TP=0;
tot=0;
[ii,jj]=find(dag==1);
for i=1:length(ii)
   
      if dag(jj(i),ii(i))==1
       if ii(i)<jj(i)
         if G(ii(i),jj(i))==1 |G(jj(i),ii(i))==1 | G(ii(i),jj(i))==2
           TP=TP+1;
         end
         tot=tot+1;
       end
      else
        tot=tot+1;
        if G(ii(i),jj(i))==1 |G(jj(i),ii(i))==1 | G(ii(i),jj(i))==2
           TP=TP+1;
        end
      end 
    
end
recall=TP./tot;


[ii,jj]=find(G==1);
n1=length(ii);
[ii,jj]=find(G==2);
n2=length(ii)/2;
precision=TP./(n1+n2);

TN=0;
[ii,jj]=find(dag==0);
for i=1:length(ii)
    if dag(jj(i),ii(i))==0
      if ii(i)<jj(i)
        if G(ii(i),jj(i))==0 & G(jj(i),ii(i))==0
          TN=TN+1;
        end
      end
    end
end
acc=(TP+TN)/(size(dag,1)*(size(dag,1)-1)/2);
        