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

function x = gen_profiles(A, n_profiles, V, theta, h, lambda, n_time_points, knockout)
% attaining artificial data for a kinetic network 
% built from a given connectivity matrix
% simulate a reaction kinetics network
if nargin < 6
	error('6 arguments required');
elseif nargin < 8
	knockout = true;
end
if nargin == 6
	n_time_points = 1; % steady state experiments
end
n = size(A, 1);
if ~isequal(size(A), [n n])
	error('A must be a square matrix');
end
if ~isscalar(n_profiles)
	error('n_profiles must be a scalar');
end
if ~isvector(V) || length(V) ~= n
	error('V must be a vector of length n')
end
if ~isequal(size(theta), [n n])
	error('theta must have the same dimensions of A')
end
if ~isequal(size(h), [n n])
	error('h must have the same dimensions of A')
end
if ~isvector(lambda) || length(lambda) ~= n
	error('lambda must be a vector of length n')
end
if ~isscalar(n_time_points) || n_time_points < 1
	error('n_time_points must be a scalar >= 1')
end

% integration routine for kinetics network
final_time = 450;
time_points = round(0:final_time/(n_time_points):final_time); % time points

n_exp = ceil(n_profiles / n_time_points);
x = zeros(n_exp * n_time_points, n);
for i = 1:n_exp
	x0 = 10 * rand(n, 1); % initial condition
    %x0=10*ones(n,1);
	if knockout
		Vcn = V;
		Vcn(mod(i-1, n) + 1) = 0;
		[tpart, xpart] = ode45(@(t, x) netw_kinet(t, x, A, Vcn, theta, h, lambda), time_points, x0);
	else
		[tpart, xpart] = ode45(@(t, x) netw_kinet(t, x, A, V, theta, h, lambda), time_points, x0);
	end
	if n_time_points > 1
		x((i-1)*n_time_points+1:i*n_time_points, :)  = xpart(2:end, :);
	else
		x(i, :) = xpart(end, :);
	end
end
if n_time_points > 1
	% return no more than n_profiles rows
	x = x(1:n_profiles, :);
end
