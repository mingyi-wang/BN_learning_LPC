function [log_Score_local] = Gauss_Score_complete_local(index,index_parents,data,T_0,T_m,v,alpha)


[n, m] = size(data); 

relevant = sort([index;index_parents]);

[log_Score_i_nom]   = Gauss_Score_complete(length(relevant),      m, v, alpha, T_0(relevant,relevant)          , T_m(relevant,relevant));

[log_Score_i_denom] = Gauss_Score_complete(length(index_parents), m, v, alpha, T_0(index_parents,index_parents), T_m(index_parents,index_parents));
    
log_Score_local = log_Score_i_nom -log_Score_i_denom;     

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [log_Score_i] = Gauss_Score_complete(n_new, m, v, alpha, T_0, T_m)

if n_new==0
log_Score_i = 0;
else
   
sum_1 = (-1)*(n_new*m)/2 * log(2*pi) + (n_new/2) * log(v/(v+m));

sum_3 = (alpha/2) * log(det(T_0)) + (-1)*(alpha+m)/2 * log(det(T_m));

log_value_a = ( alpha    *n_new/2)*log(2) + (n_new*(n_new-1)/4)*log(pi);

log_value_b = ( (alpha+m)*n_new/2)*log(2) + (n_new*(n_new-1)/4)*log(pi);
 
for i=1:n_new
    log_value_a = log_value_a + gammaln(( alpha   +1-i)/2);  % this is: log( c(n,alpha)   )
    log_value_b = log_value_b + gammaln(((alpha+m)+1-i)/2);  % this is: log( c(n,alpha+m) )
end
   
log_Score_i = sum_1 + (-1)*(log_value_a - log_value_b) + sum_3; 
 
end

return 
    
