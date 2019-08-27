function DGERocCurveNoDiag2(TrajectoryMean,idxtm,LimitX,FigNumber)

%%% TrajectoryMean --> This is the mean of the sampled DAGs during the
%%% MCMC simulation. Use function media.m to get this from a complete
%%% trajectory.
%%%
%%% idxtm --> Index of true matrix. Should be a 2 colum matrix with each line
%%% representing and edge A--->B.
%%% For example the matrix [1 2; 3 1; 4 2] means the follwing edges
%%% 1-->2; 3-->1; 4-->2
%%%
%%% LimitX is the upper range of X-axis to be presented in the graph. Note
%%% that the AUROC area under ROC curve will also be calculated with this
%%% limit. This means that you can have partial AUROC values, as AUROC_x,
%%% for any x, 0<x<=1.



%%% FigNumber is the number of the figure that will be the output.
if LimitX<=0 | LimitX>1
    display('Why?')
    why
    error('LimitX must lie between 0<LimitX<=1, type help CreateRocCurve for explanation.')
end

data_abs=TrajectoryMean;

% function to calculate where are the true edges.
true_m=DGETrueMatrix(size(data_abs),idxtm);

% function to calculate wheter there are not edges, false positives.
false_m=DGEFalseMatrixNoDiag(size(data_abs),idxtm);

% calculate the number of true edges
tot_true=sum(sum(true_m,2),1);
tot_false=sum(sum(false_m,2),1);


cnt=0;
Y=reshape(data_abs,[prod(size(data_abs)) 1]);
Y=sort(Y);
X(1,1)=Y(1,1);
novo=1;
for kk=1:size(Y,1)-1
    if Y(kk+1,1)~=Y(kk,1)
        novo=novo+1;
        X(novo,1)=Y(kk+1,1);
    end
end
         
for Step=1:size(X,1)
    cnt=cnt+1;
    found_mat=zeros(size(data_abs));
    for k=1:size(data_abs,1)
        for j=1:size(data_abs,2)
            if data_abs(k,j)>=X(Step);
                found_mat(k,j)=1;
            end
        end
    end
    %found_mat
    p_nro_corrects=true_m.*found_mat;% wheter correct coincide with found
    nro_corrects=sum(sum(p_nro_corrects,2),1);%number of correct found
    p_correct(cnt)=nro_corrects/tot_true;%relative number of corrects.
    
    p_nro_false=false_m.*found_mat;% wheter false coincide with found
    nro_false=sum(sum(p_nro_false,2),1);%number of correct found
    p_false(cnt)=nro_false/tot_false;%relative number of corrects.
end
 

p_false(end+1)=0;
p_correct(end+1)=0;
DGEPlotRocCurve2(p_false,p_correct,LimitX,FigNumber)