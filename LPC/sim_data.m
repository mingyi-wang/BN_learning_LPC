
n_rec=1;
n_genes = 1000;
edge_num=50000;
% type = 'scale-free'; % 'scale-free', 'random' or 'small-world'
%type='random';
%type='small-world';
n_time_points = 1; % number of time points for each experiment, if 1 we have steady state experiments
%n_time_points=10;
knockout = true; % whether to generate data using knockout
%knockout=false;
noise_stddev = 0.01; % standard deviation of the normal distribution of the noise added to the data matrix
max_alpha = 0.005; % maximum alpha for the calculation of pAUC(ROC)
% % FP_limit = round(n_genes / 5); % number of FP for TP_for_fixed_FP
% k = 2; % for MI_matrix and minMI_cond1var
save_all = true; % whether to save matrix, results and figures
% skip_lengthy_algorithms = n_genes > 150; % whether to skip lenghty reconstruction algorithms
script_mode = false; % whether to close Matlab at the end

% Artificial network parameters
% avg_degree = 50; % the average node degree of the artificial network
% %avg_degree=2.5;
% smallworld_p = 0.2; % the probability of link rewiring for the small-world network generation algorithm

% Other parameters
%n_measurements = round(n_genes ./ [10 5 2 1 1/2 1/5 1/10]); % vector of number of measurements. The number of experiments will be n_measurements / n_time_points
%n_measurements = round(n_genes ./ [3 1 1/2 1/5 1/30]);
n_measurements = [10 20 50 100 200 500 1000];
% m_MI = max(round(sqrt(n_measurements ./ 10)), (k+1) * ones(size(n_measurements))); % for MI_matrix; formula derived from Daub, C. O. "Analysis of integrated transcriptomics and metabolomics data - a systems biology approach" PhD thesis
% m_minMI = max(round(sqrt(n_measurements ./ 13)), (k+1) * ones(size(n_measurements))); % for minMI_cond1var
output_directory = ['../data/data_' num2str(n_genes) 'genes'];
filenames_start = [num2str(n_genes) 'g_' num2str(n_rec) 'rec_' ];
%filenames_start = [num2str(n_genes) 'g_' num2str(n_rec) 'rec_' type '_timecourse'];

% Here starts the code
if output_directory & output_directory(end) ~= '/'
	output_directory = [output_directory '/'];
end
if n_genes <= 500
	ggm_protect = 0.1;
else
	ggm_protect = 0.001;
end
n_n_measurements = numel(n_measurements);
% time_corr = zeros(n_rec, n_n_measurements);
% time_o1pcorr = zeros(n_rec, n_n_measurements);
% time_o2pcorr = zeros(n_rec, n_n_measurements);
% time_ggm = zeros(n_rec, n_n_measurements);
% time_MI = zeros(n_rec, n_n_measurements);
% time_minMI = zeros(n_rec, n_n_measurements);
% time_DPI = zeros(n_rec, n_n_measurements);
% AUC_ROC_corr = zeros(n_rec, n_n_measurements);
% AUC_ROC_o1pcorr = zeros(n_rec, n_n_measurements);
% AUC_ROC_o2pcorr = zeros(n_rec, n_n_measurements);
% AUC_ROC_ggm = zeros(n_rec, n_n_measurements);
% AUC_ROC_MI = zeros(n_rec, n_n_measurements);
% AUC_ROC_minMI = zeros(n_rec, n_n_measurements);
% AUC_ROC_DPI = zeros(n_rec, n_n_measurements);
% pAUC_ROC_corr = zeros(n_rec, n_n_measurements);
% pAUC_ROC_o1pcorr = zeros(n_rec, n_n_measurements);
% pAUC_ROC_o2pcorr = zeros(n_rec, n_n_measurements);
% pAUC_ROC_ggm = zeros(n_rec, n_n_measurements);
% pAUC_ROC_MI = zeros(n_rec, n_n_measurements);
% pAUC_ROC_minMI = zeros(n_rec, n_n_measurements);
% pAUC_ROC_DPI = zeros(n_rec, n_n_measurements);
% AUC_PvsR_corr = zeros(n_rec, n_n_measurements);
% AUC_PvsR_o1pcorr = zeros(n_rec, n_n_measurements);
% AUC_PvsR_o2pcorr = zeros(n_rec, n_n_measurements);
% AUC_PvsR_ggm = zeros(n_rec, n_n_measurements);
% AUC_PvsR_MI = zeros(n_rec, n_n_measurements);
% AUC_PvsR_minMI = zeros(n_rec, n_n_measurements);
% AUC_PvsR_DPI = zeros(n_rec, n_n_measurements);
% TP_corr = zeros(n_rec, n_n_measurements);
% TP_o1pcorr = zeros(n_rec, n_n_measurements);
% TP_o2pcorr = zeros(n_rec, n_n_measurements);
% TP_ggm = zeros(n_rec, n_n_measurements);
% TP_MI = zeros(n_rec, n_n_measurements);
% TP_minMI = zeros(n_rec, n_n_measurements);
% TP_DPI = zeros(n_rec, n_n_measurements);

% initializeR({'Rmatlab' '--no-save' '--no-restore'})
% callR('library', 'GeneNet')

for i = 1:n_rec
	i
% 	switch type
% 		case 'scale-free'
% 			%A = scalefree_BA(n_genes, avg_degree);
            edges=0;
            while edges<edge_num
              [A, ord,edges]=mk_rnd_dag(n_genes,100,edge_num);               
            end
% 		case 'random'
% 			A = random_net(n_genes, avg_degree);
% 		case 'small-world'
% 			A = smallworld_net(n_genes, avg_degree, smallworld_p);
% 	end
	%A_undir = A | A'; % the undirected matrix used to measure TP, FP,...

	% add random sign to each edge (-1 or 1)
	A = sign(A .* (rand(n_genes) - 0.5));
	% choose the basal values of transcription of each gene
	%V = rand(n_genes, 1);
    V=1*ones(n_genes,1);
	% choose the half-life coefficients
	%theta = 10 * rand(n_genes);
    theta=100*ones(n_genes);
	% choose Hill coefficients in 1:4
    % h = ceil(4 * rand(n_genes));
    %h= 1.5*ones(n_genes);
    h= 1*ones(n_genes);
	% choose the degradation rates
	lambda = 0.1 * rand(n_genes, 1);

	if save_all
		%save([output_directory filenames_start '_A_' num2str(i) '.mat'], 'A', 'V', 'theta', 'h', 'lambda')
        dlmwrite([output_directory filenames_start 'A_' num2str(i) '.txt'], A,'\t');
	end

	for j = 2:n_n_measurements
		j
% 		if n_rec > 1
% 			plotCurves = 0;
% 		else
% 			plotCurves(1) = figure;
% 			hold all
% 			plotCurves(2) = figure;
% 			hold all
% 		end

		x = gen_profiles(A, n_measurements(j), V, theta, h, lambda, n_time_points, knockout);
		% Add Gaussian noise to x
		x = x + noise_stddev * randn(size(x));
		x(x < 0) = 0;
		if save_all
			%save([output_directory filenames_start '_x_' num2str(i) '_' num2str(n_measurements(j)) 'measurements.mat'], 'x')
            dlmwrite([output_directory filenames_start 'x_' num2str(i) '_' num2str(n_measurements(j)) 'measurements.txt'], x,'\t');
		end
		if n_time_points > 1
			% rank transformation of the columns of x
			xr = zeros(size(x));
			for u = 1:n_genes
				xr(:, u) = rank_transf(x(:, u));
			end
			xr = xr / size(x, 1);
			clear u
		else
			xr = x;
		end

% 		tic
% 		W_corr = abs(corr(x));
% 		time_corr(i, j) = toc;
% 		if save_all
% 			save([output_directory filenames_start '_W_corr_' num2str(i) '_' num2str(n_measurements(j)) 'measurements.mat'], 'W_corr')
% 		end
% 		[AUC_ROC_corr(i, j) pAUC_ROC_corr(i, j) AUC_PvsR_corr(i,j)] = AUCs(W_corr, A_undir, plotCurves, max_alpha);
% 		TP_corr(i, j) = TP_for_fixed_FP(W_corr, A_undir, FP_limit);
% 
% 		tic
% 		W_o1pcorr = abs(min_o1pcorr(x));
% 		time_o1pcorr(i, j) = toc;
% 		if save_all
% 			save([output_directory filenames_start '_W_o1pcorr_' num2str(i) '_' num2str(n_measurements(j)) 'measurements.mat'], 'W_o1pcorr')
% 		end
% 		[AUC_ROC_o1pcorr(i, j) pAUC_ROC_o1pcorr(i, j) AUC_PvsR_o1pcorr(i, j)] = AUCs(W_o1pcorr, A_undir, plotCurves, max_alpha);
% 		TP_o1pcorr(i, j) = TP_for_fixed_FP(W_o1pcorr, A_undir, FP_limit);
% 
% 		if ~skip_lengthy_algorithms
% 			tic
% 			W_o2pcorr = abs(min_o2pcorr(x));
% 			time_o2pcorr(i, j) = toc;
% 			if save_all
% 				save([output_directory filenames_start '_W_o2pcorr_' num2str(i) '_' num2str(n_measurements(j)) 'measurements.mat'], 'W_o2pcorr')
% 			end
% 			[AUC_ROC_o2pcorr(i, j) pAUC_ROC_o2pcorr(i, j) AUC_PvsR_o2pcorr(i, j)] = AUCs(W_o2pcorr, A_undir, plotCurves, max_alpha);
% 			TP_o2pcorr(i, j) = TP_for_fixed_FP(W_o2pcorr, A_undir, FP_limit);
% 		end
% 
% 		tic
% 		W_ggm = abs(callNamedR('ggm.estimate.pcor', x, {'protect', ggm_protect}));
% 		time_ggm(i, j) = toc;
% 		if nnz(isnan(W_ggm))
% 			warning('W_ggm contains NaNs.')
% 		end
% 		if nnz(W_ggm > 1)
% 			warning('W_ggm contains elements > 1.')
% 		end
% 		if save_all
% 			save([output_directory filenames_start '_W_ggm_' num2str(i) '_' num2str(n_measurements(j)) 'measurements.mat'], 'W_ggm')
% 		end
% 		[AUC_ROC_ggm(i, j) pAUC_ROC_ggm(i, j) AUC_PvsR_ggm(i, j)] = AUCs(W_ggm, A_undir, plotCurves, max_alpha);
% 		TP_ggm(i, j) = TP_for_fixed_FP(W_ggm, A_undir, FP_limit);
% 
% 		tic
% 		W_MI = MI_matrix(xr', m_MI(j), k);
% 		time_MI(i, j) = toc;
% 		if save_all
% 			save([output_directory filenames_start '_W_MI_' num2str(i) '_' num2str(n_measurements(j)) 'measurements.mat'], 'W_MI')
% 		end
% 		[AUC_ROC_MI(i, j) pAUC_ROC_MI(i, j) AUC_PvsR_MI(i ,j)] = AUCs(W_MI, A_undir, plotCurves, max_alpha);
% 		TP_MI(i, j) = TP_for_fixed_FP(W_MI, A_undir, FP_limit);
% 
% 		if ~skip_lengthy_algorithms
% 			tic
% 			W_minMI = minMI_cond1var(xr', m_minMI(j), k);
% 			time_minMI(i, j) = toc;
% 			if save_all
% 				save([output_directory filenames_start '_W_minMI_' num2str(i) '_' num2str(n_measurements(j)) 'measurements.mat'], 'W_minMI')
% 			end
% 			[AUC_ROC_minMI(i, j) pAUC_ROC_minMI(i, j) AUC_PvsR_minMI(i, j)] = AUCs(W_minMI, A_undir, plotCurves, max_alpha);
% 			TP_minMI(i, j) = TP_for_fixed_FP(W_minMI, A_undir, FP_limit);
% 		end
% 
% 		tic
% 		W_DPI = DPI_pruning(W_MI, 0, 0.1);
% 		time_DPI(i, j) = toc + time_MI(i, j);
% 		if save_all
% 			save([output_directory filenames_start '_W_DPI_' num2str(i) '_' num2str(n_measurements(j)) 'measurements.mat'], 'W_DPI')
% 		end
% 		[AUC_ROC_DPI(i, j) pAUC_ROC_DPI(i, j) AUC_PvsR_DPI(i, j)] = AUCs(W_DPI, A_undir, plotCurves, max_alpha);
% 		TP_DPI(i, j) = TP_for_fixed_FP(W_DPI, A_undir, FP_limit);
% 
% 		if plotCurves
% 			figure(plotCurves(1))
% 			title(['ROC curves for a ' type ' network of ' num2str(n_genes) ' genes and data from ' num2str(n_measurements(j)) ' measurements'])
% 			legend('R', 'R_{C1}', 'R_{C2}', 'I', 'I_C', 'I_{DPI}', 'R_{Call}', 'Location', 'SouthEast')
% 			if save_all
% 				saveas(gcf, [output_directory filenames_start '_ROC_' num2str(n_measurements(j)) 'measurements.fig'])
% 			end
% 
% 			figure(plotCurves(2))
% 			title(['PvsR curves for a ' type ' network of ' num2str(n_genes) ' genes and data from ' num2str(n_measurements(j)) ' measurements'])
% 			legend('R', 'R_{C1}', 'R_{C2}', 'I', 'I_C', 'I_{DPI}', 'R_{Call}', 'Location', 'SouthEast')
% 			if save_all
% 				saveas(gcf, [output_directory filenames_start '_PvsR_' num2str(n_measurements(j)) 'measurements.fig'])
% 			end
% 
% 			figure
% 			plot(W_MI, W_corr, '.k', 'LineStyle', 'none')
% 			xlabel('Mutual informations')
% 			ylabel('Pearson correlations')
% 			hold all
% 			Y = [0:0.05:0.95];
% 			plot(-log(sqrt(ones(size(Y)) - Y .^ 2)), Y, 'LineWidth', 2)
% 			if save_all
% 				saveas(gcf, [output_directory filenames_start '_MIvscorr_' num2str(n_measurements(j)) 'measurements.fig'])
% 			end
% 		end
	end
end

% clear A A_undir V theta h lambda x xr W_*
% 
% time_corr_mean = mean(time_corr, 1);
% time_corr_std = std(time_corr, 0, 1);
% time_o1pcorr_mean = mean(time_o1pcorr, 1);
% time_o1pcorr_std = std(time_o1pcorr, 0, 1);
% time_o2pcorr_mean = mean(time_o2pcorr, 1);
% time_o2pcorr_std = std(time_o2pcorr, 0, 1);
% time_ggm_mean = mean(time_ggm, 1);
% time_ggm_std = std(time_ggm, 0, 1);
% time_MI_mean = mean(time_MI, 1);
% time_MI_std = std(time_MI, 0, 1);
% time_minMI_mean = mean(time_minMI, 1);
% time_minMI_std = std(time_minMI, 0, 1);
% time_DPI_mean = mean(time_DPI, 1);
% time_DPI_std = std(time_DPI, 0, 1);
% 
% AUC_ROC_corr_mean = mean(AUC_ROC_corr, 1);
% AUC_ROC_corr_std = std(AUC_ROC_corr, 0, 1);
% AUC_ROC_o1pcorr_mean = mean(AUC_ROC_o1pcorr, 1);
% AUC_ROC_o1pcorr_std = std(AUC_ROC_o1pcorr, 0, 1);
% AUC_ROC_o2pcorr_mean = mean(AUC_ROC_o2pcorr, 1);
% AUC_ROC_o2pcorr_std = std(AUC_ROC_o2pcorr, 0, 1);
% AUC_ROC_ggm_mean = mean(AUC_ROC_ggm, 1);
% AUC_ROC_ggm_std = std(AUC_ROC_ggm, 0, 1);
% AUC_ROC_MI_mean = mean(AUC_ROC_MI, 1);
% AUC_ROC_MI_std = std(AUC_ROC_MI, 0, 1);
% AUC_ROC_minMI_mean = mean(AUC_ROC_minMI, 1);
% AUC_ROC_minMI_std = std(AUC_ROC_minMI, 0, 1);
% AUC_ROC_DPI_mean = mean(AUC_ROC_DPI, 1);
% AUC_ROC_DPI_std = std(AUC_ROC_DPI, 0, 1);
% 
% pAUC_ROC_corr_mean = mean(pAUC_ROC_corr, 1);
% pAUC_ROC_corr_std = std(pAUC_ROC_corr, 0, 1);
% pAUC_ROC_o1pcorr_mean = mean(pAUC_ROC_o1pcorr, 1);
% pAUC_ROC_o1pcorr_std = std(pAUC_ROC_o1pcorr, 0, 1);
% pAUC_ROC_o2pcorr_mean = mean(pAUC_ROC_o2pcorr, 1);
% pAUC_ROC_o2pcorr_std = std(pAUC_ROC_o2pcorr, 0, 1);
% pAUC_ROC_ggm_mean = mean(pAUC_ROC_ggm, 1);
% pAUC_ROC_ggm_std = std(pAUC_ROC_ggm, 0, 1);
% pAUC_ROC_MI_mean = mean(pAUC_ROC_MI, 1);
% pAUC_ROC_MI_std = std(pAUC_ROC_MI, 0, 1);
% pAUC_ROC_minMI_mean = mean(pAUC_ROC_minMI, 1);
% pAUC_ROC_minMI_std = std(pAUC_ROC_minMI, 0, 1);
% pAUC_ROC_DPI_mean = mean(pAUC_ROC_DPI, 1);
% pAUC_ROC_DPI_std = std(pAUC_ROC_DPI, 0, 1);
% 
% AUC_PvsR_corr_mean = mean(AUC_PvsR_corr, 1);
% AUC_PvsR_corr_std = std(AUC_PvsR_corr, 0, 1);
% AUC_PvsR_o1pcorr_mean = mean(AUC_PvsR_o1pcorr, 1);
% AUC_PvsR_o1pcorr_std = std(AUC_PvsR_o1pcorr, 0, 1);
% AUC_PvsR_o2pcorr_mean = mean(AUC_PvsR_o2pcorr, 1);
% AUC_PvsR_o2pcorr_std = std(AUC_PvsR_o2pcorr, 0, 1);
% AUC_PvsR_ggm_mean = mean(AUC_PvsR_ggm, 1);
% AUC_PvsR_ggm_std = std(AUC_PvsR_ggm, 0, 1);
% AUC_PvsR_MI_mean = mean(AUC_PvsR_MI, 1);
% AUC_PvsR_MI_std = std(AUC_PvsR_MI, 0, 1);
% AUC_PvsR_minMI_mean = mean(AUC_PvsR_minMI, 1);
% AUC_PvsR_minMI_std = std(AUC_PvsR_minMI, 0, 1);
% AUC_PvsR_DPI_mean = mean(AUC_PvsR_DPI, 1);
% AUC_PvsR_DPI_std = std(AUC_PvsR_DPI, 0, 1);
% 
% TP_corr_mean = mean(TP_corr, 1);
% TP_corr_std = std(TP_corr, 0, 1);
% TP_o1pcorr_mean = mean(TP_o1pcorr, 1);
% TP_o1pcorr_std = std(TP_o1pcorr, 0, 1);
% TP_o2pcorr_mean = mean(TP_o2pcorr, 1);
% TP_o2pcorr_std = std(TP_o2pcorr, 0, 1);
% TP_ggm_mean = mean(TP_ggm, 1);
% TP_ggm_std = std(TP_ggm, 0, 1);
% TP_MI_mean = mean(TP_MI, 1);
% TP_MI_std = std(TP_MI, 0, 1);
% TP_minMI_mean = mean(TP_minMI, 1);
% TP_minMI_std = std(TP_minMI, 0, 1);
% TP_DPI_mean = mean(TP_DPI, 1);
% TP_DPI_std = std(TP_DPI, 0, 1);
% 
% if n_n_measurements > 1
% 	figure
% 	hold all
% 	plot(n_measurements, time_corr_mean, '-db')
% 	plot(n_measurements, time_o1pcorr_mean, '-v', 'Color', [0 0.5 0])
% 	if ~skip_lengthy_algorithms
% 		plot(n_measurements, time_o2pcorr_mean, '-^r')
% 	end
% 	plot(n_measurements, time_ggm_mean, '-p', 'Color', [0 0.75 0.75])
% 	plot(n_measurements, time_MI_mean, '-s', 'Color', [0.75 0 0.75])
% 	if ~skip_lengthy_algorithms
% 		plot(n_measurements, time_minMI_mean, '-o', 'Color', [0.75 0.75 0])
% 	end
% 	plot(n_measurements, time_DPI_mean, '-*', 'Color', [0.25 0.25 0.25])
% 	xlabel('N. of measurements')
% 	ylabel('Elapsed time (s)')
% 	title(['Elapsed time for ' type ' networks of ' num2str(n_genes) ' genes'])
% 	if ~skip_lengthy_algorithms
% 		legend('R', 'R_{C1}', 'R_{C2}', 'R_{Call}', 'I', 'I_C', 'I_{DPI}', 'Location', 'West')
% 	else
% 		legend('R', 'R_{C1}', 'R_{Call}', 'I', 'I_{DPI}', 'Location', 'West')
% 	end
% 	if save_all
% 		saveas(gcf, [output_directory filenames_start '_time.fig'])
% 	end
% 
% 	figure
% 	hold all
% 	plot(n_measurements, AUC_ROC_corr_mean, '-db')
% 	plot(n_measurements, AUC_ROC_o1pcorr_mean, '-v', 'Color', [0 0.5 0])
% 	if ~skip_lengthy_algorithms
% 		plot(n_measurements, AUC_ROC_o2pcorr_mean, '-^r')
% 	end
% 	plot(n_measurements, AUC_ROC_ggm_mean, '-p', 'Color', [0 0.75 0.75])
% 	plot(n_measurements, AUC_ROC_MI_mean, '-s', 'Color', [0.75 0 0.75])
% 	if ~skip_lengthy_algorithms
% 		plot(n_measurements, AUC_ROC_minMI_mean, '-o', 'Color', [0.75 0.75 0])
% 	end
% 	plot(n_measurements, AUC_ROC_DPI_mean, '-*', 'Color', [0.25 0.25 0.25])
% 	xlabel('N. of measurements')
% 	ylabel('AUC(ROC)')
% 	ylim([0 1])
% 	title(['AUC(ROC) for ' type ' networks of ' num2str(n_genes) ' genes'])
% 	if ~skip_lengthy_algorithms
% 		legend('R', 'R_{C1}', 'R_{C2}', 'R_{Call}', 'I', 'I_C', 'I_{DPI}', 'Location', 'SouthEast')
% 	else
% 		legend('R', 'R_{C1}', 'R_{Call}', 'I', 'I_{DPI}', 'Location', 'SouthEast')
% 	end
% 	if save_all
% 		saveas(gcf, [output_directory filenames_start '_AUC_ROC.fig'])
% 	end
% 
% 	figure
% 	hold all
% 	plot(n_measurements, pAUC_ROC_corr_mean, '-db')
% 	plot(n_measurements, pAUC_ROC_o1pcorr_mean, '-v', 'Color', [0 0.5 0])
% 	if ~skip_lengthy_algorithms
% 		plot(n_measurements, pAUC_ROC_o2pcorr_mean, '-^r')
% 	end
% 	plot(n_measurements, pAUC_ROC_ggm_mean, '-p', 'Color', [0 0.75 0.75])
% 	plot(n_measurements, pAUC_ROC_MI_mean, '-s', 'Color', [0.75 0 0.75])
% 	if ~skip_lengthy_algorithms
% 		plot(n_measurements, pAUC_ROC_minMI_mean, '-o', 'Color', [0.75 0.75 0])
% 	end
% 	plot(n_measurements, pAUC_ROC_DPI_mean, '-*', 'Color', [0.25 0.25 0.25])
% 	xlabel('N. of measurements')
% 	ylabel('pAUC(ROC)')
% 	title(['pAUC(ROC) (alpha<=' num2str(max_alpha) ') for ' type ' networks of ' num2str(n_genes) ' genes'])
% 	if ~skip_lengthy_algorithms
% 		legend('R', 'R_{C1}', 'R_{C2}', 'R_{Call}', 'I', 'I_C', 'I_{DPI}', 'Location', 'SouthEast')
% 	else
% 		legend('R', 'R_{C1}', 'R_{Call}', 'I', 'I_{DPI}', 'Location', 'SouthEast')
% 	end
% 	if save_all
% 		saveas(gcf, [output_directory filenames_start '_pAUC_ROC.fig'])
% 	end
% 
% 	figure
% 	hold all
% 	plot(n_measurements, AUC_PvsR_corr_mean, '-db')
% 	plot(n_measurements, AUC_PvsR_o1pcorr_mean, '-v', 'Color', [0 0.5 0])
% 	if ~skip_lengthy_algorithms
% 		plot(n_measurements, AUC_PvsR_o2pcorr_mean, '-^r')
% 	end
% 	plot(n_measurements, AUC_PvsR_ggm_mean, '-p', 'Color', [0 0.75 0.75])
% 	plot(n_measurements, AUC_PvsR_MI_mean, '-s', 'Color', [0.75 0 0.75])
% 	if ~skip_lengthy_algorithms
% 		plot(n_measurements, AUC_PvsR_minMI_mean, '-o', 'Color', [0.75 0.75 0])
% 	end
% 	plot(n_measurements, AUC_PvsR_DPI_mean, '-*', 'Color', [0.25 0.25 0.25])
% 	xlabel('N. of measurements')
% 	ylabel('AUC(PvsR)')
% 	ylim([0 0.7])
% 	title(['AUC(PvsR) for ' type ' networks of ' num2str(n_genes) ' genes'])
% 	if ~skip_lengthy_algorithms
% 		legend('R', 'R_{C1}', 'R_{C2}', 'R_{Call}', 'I', 'I_C', 'I_{DPI}', 'Location', 'NorthWest')
% 	else
% 		legend('R', 'R_{C1}', 'R_{Call}', 'I', 'I_{DPI}', 'Location', 'NorthWest')
% 	end
% 	if save_all
% 		saveas(gcf, [output_directory filenames_start '_AUC_PvsR.fig'])
% 	end
% 
% 	figure
% 	hold all
% 	plot(n_measurements, TP_corr_mean, '-db')
% 	plot(n_measurements, TP_o1pcorr_mean, '-v', 'Color', [0 0.5 0])
% 	if ~skip_lengthy_algorithms
% 		plot(n_measurements, TP_o2pcorr_mean, '-^r')
% 	end
% 	plot(n_measurements, TP_ggm_mean, '-p', 'Color', [0 0.75 0.75])
% 	plot(n_measurements, TP_MI_mean, '-s', 'Color', [0.75 0 0.75])
% 	if ~skip_lengthy_algorithms
% 		plot(n_measurements, TP_minMI_mean, '-o', 'Color', [0.75 0.75 0])
% 	end
% 	plot(n_measurements, TP_DPI_mean, '-*', 'Color', [0.25 0.25 0.25])
% 	xlabel('N. of measurements')
% 	ylabel('TP')
% 	title(['TP for FP=' num2str(FP_limit) ' for ' type ' networks of ' num2str(n_genes) ' genes'])
% 	if ~skip_lengthy_algorithms
% 		legend('R', 'R_{C1}', 'R_{C2}', 'R_{Call}', 'I', 'I_C', 'I_{DPI}', 'Location', 'SouthEast')
% 	else
% 		legend('R', 'R_{C1}', 'R_{Call}', 'I', 'I_{DPI}', 'Location', 'SouthEast')
% 	end
% 	if save_all
% 		saveas(gcf, [output_directory filenames_start '_TP.fig'])
% 	end
% end
% 
% if save_all
% 	save([output_directory filenames_start '.mat'])
% end
% if script_mode
% 	quit
% end
