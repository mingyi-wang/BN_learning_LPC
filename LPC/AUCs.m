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

% Compute the AUC(ROC) (Area Under the ROC Curve), pAUC(ROC) (partial AUC(ROC))
% and AUC(PvsR) (Area Under PvsR Curve) and plot the ROC (Receiver Operating 
% Characteristic) and PvsR (Precision vs Recall) curve for an undirected 
% network prediction.
%
% W             square matrix of edge weights of the predicted network
% A             square boolean matrix being the incidence matrix for the real network
% plotCurves    0 to not plot curves, vector of two figure handles otherwise
% max_alpha     maximum alpha for the calculation of pAUC
function [AUC_ROC pAUC_ROC AUC_PvsR] = AUCs(W, A, plotCurves, max_alpha)
if nargin < 2
    error('2 arguments required');
elseif nargin < 4
    max_alpha = 0.05;
end
if nargin == 2
    plotCurves = 0;
end
if ndims(W) > 2
    error('W must be a matrix');
end
n = size(W, 1);
if n ~= size(W, 2)
    error('W must be a square matrix');
end
if ~isequal(size(W), size(A))
    error('A must have the same dimensions of W');
end

% Set to 0 the diagonal elements of W
for i = 1:n
    W(i, i) = 0;
end
% Set to 0 the NaN in W
W(isnan(W)) = 0;

% The possible threshold are the unique values of W
thresholds = sort(unique(W(:)), 'descend');
% Limit the number of thresholds to save time
max_n_thresholds = 10000;
if numel(thresholds) > max_n_thresholds
    thresholds = thresholds(round(1:(numel(thresholds) - 1) / (max_n_thresholds - 1):numel(thresholds)));
end

% # of real positives and negatives (if A is simmetric, i.e. the real graph is
% undirected, then all the following counts will be doubled)
RealPos = nnz(A);
n_edges = n^2 - n;
RealNeg = n_edges - RealPos;

% Calculate TP and FP for all the possible thresholds
TP = zeros(1, numel(thresholds));
FP = zeros(1, numel(thresholds));
for i = 1:numel(thresholds)
    Arec = W >= thresholds(i);
    % # of true positives
    TP(i) = nnz(A & Arec);
    % # of false positives = # predicted positive - TP
    for j=1:length(Arec)
        Arec(j,j)=0;
    end 
    FP(i) = nnz(Arec) - TP(i);
    
end

sensitivity = TP ./ RealPos;
% alpha = 1 - specificity
alpha = FP ./ RealNeg;
precision = TP ./ (TP + FP);

max_index = find(alpha <= max_alpha, 1, 'last');

if plotCurves
    % plot the ROC curve
    figure(plotCurves(1))
    
    plot(alpha, sensitivity);
    axis equal
    axis([0 1 0 1])
    xlabel('1 - specificity');
    ylabel('sensitivity');
    title('ROC curve');

    % plot the PvsR curve
    figure(plotCurves(2))
    
    plot(sensitivity, precision);
    axis equal
    xlabel('recall');
    ylabel('precision');
    title('precision vs recall');
end

if max_index > 1
    pAUC_ROC = trapz(alpha(1:max_index), sensitivity(1:max_index));
else
    pAUC_ROC = 0;
end
if numel(sensitivity) > 1
    % Calculate the Area Under the ROC Curve
    AUC_ROC = trapz(alpha, sensitivity);
    % Calculate the Area Under PvsR Curve
    AUC_PvsR = trapz(sensitivity, precision);
else
    AUC_ROC = 0;
    AUC_PvsR = 0;
end
