% Copyright (C) 2006, 2007 Nicola Soranzo <soranzo@sissa.it>
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

% W weight matrix n * n
% A true boolean matrix
function TP = TP_for_fixed_FP(W, A, FP_limit)
if nargin < 3
    error('3 arguments required');
end
n = size(W, 1);
% Set to 0 the NaN in W
W(isnan(W)) = 0;
% sort the upper triangular part of W in descending order
[edge_sorted edge_indexes] = sort(reshape(triu(W, 1), n^2, 1), 'descend');

FP = 0;
i = 1;
while FP <= FP_limit
    if ~A(edge_indexes(i))
        FP = FP + 1;
    end
    i = i + 1;
end
TP = i - 1 - FP;
