%% example how to get the relevance networks using the pearson correlation

% in this example we get the relevance network pearson correlation for the
% 1st data set that was generated with multivariate gaussian distribution
% and with observations only.

load Gaussian_Obs
RelPe=abs(corrcoef(data_gaussian_obs{1}));

