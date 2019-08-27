% Copyright (C) 2006 Claudio Altafini <altafini@sissa.it>
% Copyright (C) 2006,2007 Nicola Soranzo <soranzo@sissa.it>
%
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the Free 
% Software Foundation; either version 2 of the License, or (at your option) 
% any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
% for more details.

% You should have received a copy of the GNU General Public License 
% along with this program; see the file COPYING. If not, write to the Free 
% Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
% 02110-1301, USA.

% Function to generate reaction kinetics according to
% a given adjacency matrix A
function dxdt = netw_kinet(t, x, A, V, theta, h, lambda)
n = length(x);
dxdt = V;
indegrees = sum(A ~= 0, 2);
for i = 1:n
    if V(i) & indegrees(i)
        for j = 1:n
            if A(i, j) % there is an arc (i, j)
                % choose Hill parameter
                if A(i, j) > 0 % activation
                     x_j_to_h_ij = x(j)^h(i, j);
                     dxdt(i) = dxdt(i) * (1 + x_j_to_h_ij / (x_j_to_h_ij + theta(i, j)^h(i, j)));
 %                   dxdt(i)=dxdt(i)*(1+x(j)/x(j)+0.01))^1.5;
                else % inhibition 
                    theta_ij_to_h_ij = theta(i, j)^h(i, j);
                    dxdt(i) = dxdt(i) * (theta_ij_to_h_ij / (x(j)^h(i, j) + theta_ij_to_h_ij));
                end
            end
        end
    end
end
dxdt = dxdt - lambda .* x; % decay rate
