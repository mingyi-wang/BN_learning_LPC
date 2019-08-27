function [score]=MCMCOmatrix2score(bayes)
bayes_sort=sort(bayes,2);
bayes_sum=sum(bayes_sort,2);
idx=find(bayes_sum==0);
bayes_sum(idx)=1;
bayes_log=log(bayes_sum);
score=sum(bayes_log);
