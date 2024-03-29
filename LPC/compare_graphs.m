function [num_a,num_a1,num_true,num_mo,num_ea,num_ma,num_wo,SHD]=compare_graphs(G,eg)
% G is the DAG which want to be compared, G(x,y)=1 means x->y; G(x,y)=2 and G(y,x)=2 means x-y
% eg is the CPAG: eg(x,y)=2 x->y  eg(x,y)=3 x-y
% 
% Structure of Hamming distance calculation
%
% num_a: the edge numbers of orignial graph
% num_a1: the edge numbers of orignial graph


num_a=size(find(eg~=0),1);
num_a1=size(find(G==2),1)/2+size(find(G==1),1);

[II JJ]=find(G==1);
num_true=0;   %True prediction
for i=1:length(II)
  if (eg(II(i),JJ(i))==2) 
    num_true=num_true+1;
  end
end

num_mo=0;     %Missing Orientation
[II JJ]=find(G==2); %this is the directed 
for i=1:length(II)
  if II(i)<JJ(i)
    if (eg(II(i),JJ(i))==2)  | (eg(JJ(i),II(i))==2)       
      num_mo=num_mo+1;       
    end
  end
end
[II JJ]=find(eg==1); %this is the directed in true eg
for i=1:length(II)
  if II(i)<JJ(i)
    if (G(II(i),JJ(i))==3)  | (G(JJ(i),II(i))==3)       
      num_mo=num_mo+1;       
    end
  end
end


%Extra Arc
[II JJ]=find(G==1|G==2);
num_ea=0;
for i=1:length(II)
  if G(II(i),JJ(i))==2  % having no direction
    if II(i)<JJ(i)
       if eg(II(i),JJ(i))==0 & eg(JJ(i),II(i))==0
          num_ea=num_ea+1;
       end
    end
  else
     if eg(II(i),JJ(i))==0 & eg(JJ(i),II(i))==0
        num_ea=num_ea+1;    
     end
  end
end

%Missing Arcs
num_ma=0;
[II JJ]=find(eg~=0);
for i=1:length(II)
    if G(II(i),JJ(i))==0 & G(JJ(i),II(i))==0
        num_ma=num_ma+1;
    end    
end


%Wrongly oriented
num_wo=0;
 num_edge=0;
[II JJ]=find(eg==2);
for i=1:length(II)
    if G(JJ(i),II(i))==1
       num_wo=num_wo+1;
    end    
    if G(JJ(i),II(i))==1 | G(II(i),JJ(i))==1
       num_edge=num_edge+1;
    end    
end
ratio=num_wo/num_edge
SHD=num_wo+num_ma+num_ea+num_mo;
%fprintf('num:%d, true:%d, M.O.:%d, E.A.:%d, M.A.:%d, W.O.:%d\n',num_a,num_true,num_mo,num_ea,num_ma,num_wo);

       