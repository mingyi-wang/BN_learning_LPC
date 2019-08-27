type='scale-free';
for i=1:10
    
    fname=sprintf('%s%s%s%d%s','H:\buf\NetRec\data\100g_10rec_',type,'_A_',i,'.txt');
    fname1=sprintf('%s%s%s%d%s','H:\buf\NetRec\data\100g_10rec_',type,'_A_',i,'.edg');
    fid=fopen(fname1,'w+');
    A=dlmread(fname,'\t');
    [ii,jj]=size(A);
    for k=1:ii
        for l=1:jj
            if A(k,l)==1 | A(k,l)==-1
                fprintf(fid,'%d\t%d\n',k,l);
            end
        end
    end
    fclose(fid);
end