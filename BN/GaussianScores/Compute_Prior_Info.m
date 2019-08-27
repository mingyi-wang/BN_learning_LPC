function [T_0, T_m, v, alpha] = Compute_Prior_Info(data)


[n, m] = size(data);%n=nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THESE VALUES HAVE TO BE DEFINED BY A USER: v, alpha, my_0, and B 
v     = 1;
alpha = n+2;
my_0  = mean(mean(data'))*ones(n,1);
v_vec = mean(var(data'))*ones(n,1); 
B     = zeros(n,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE PRECISION-MATRIX W
W =  1/v_vec(1);
for i = 1:(n-1)
    b_vec = B(1:i,(i+1));
    W = [W + (b_vec * b_vec')/v_vec(i+1) , (-1)*b_vec/v_vec(i+1) ; (-1)*b_vec'/v_vec(i+1), 1/v_vec(i+1)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE T_0 and T_m

T_0   = (v*(alpha-n-1))/(v+1) * inv(W);   

S_m = zeros(n,n);
for i=1:m    
S_m = S_m + (data(:,i) - mean(data')') * (data(:,i) - mean(data')')';
end

T_m = T_0 + S_m + (v * m)/(v+m) * (my_0 - mean(data')') * (my_0 - mean(data')')';
