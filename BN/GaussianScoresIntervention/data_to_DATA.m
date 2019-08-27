function [DATA] = data_to_DATA(data,vector)

[n_nodes, n_obs] = size(data);

for i=1:n_nodes
indicis = find(vector~=i & vector~=(-i));
DATA{i} = data(:,indicis);
end

return



