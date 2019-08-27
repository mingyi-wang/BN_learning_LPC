load ('..\data\ecoli\pdag_ecoli_udag_ord3.mat');
load('..\data\ecoli\pdag_ecoli_pdag_ord3.mat');
[ii,jj]=find(pdag~=0);


fid=fopen('..\data\ecoli\ecoli_net.txt','w');
fid1=fopen('..\data\ecoli\ecoli_genes.txt','w');
load('..\data\ecoli\compendium_E_coli_v3_Build_1.mat');
load('..\data\ecoli\tfs\tfs.mat');
load('..\data\ecoli\reg\reg.mat');
largeA=zeros(7312,7312);
for i=1:length(reg.A)
  for j=1:length(reg.A)
   if reg.A(i,j)==1
      largeA(reg.gidx_7312(i),reg.gidx_7312(j))=1;
   end
  end
end
for i=1:length(ii)
  if pdag(ii(i),jj(i))==-1
    fprintf(fid,'%s\t%s\t%d\t%d\t%f\t%d',compendium.genes{compendium.gidx(ii(i))},compendium.genes{compendium.gidx(jj(i))},ii(i),jj(i),zMin(ii(i),jj(i)),pdag(ii(i),jj(i)));
    if largeA(compendium.gidx(ii(i)),compendium.gidx(jj(i)))==1 
        fprintf(fid,'\t%d\n',-1);
    elseif largeA(compendium.gidx(jj(i)),compendium.gidx(ii(i)))==1 
        fprintf(fid,'\t%d\n',1);
    else
        fprintf(fid,'\t%d\n',0);
    end 
  end
  
   i
end
for i=1:length(compendium.gidx)
   
   fprintf(fid1,'%d\t%s\t',i,compendium.genes{compendium.gidx(i)});
   if ismember(compendium.gidx(i),tfs.gidx_7312) 
      fprintf(fid1,'%d',1); 
   else
      fprintf(fid1,'%d',0);  
   end
   degree=length(find(undirected_G(:,i)~=0));
   fprintf(fid1,'\t%d\n',degree);
end
fclose(fid);
fclose(fid1);
    