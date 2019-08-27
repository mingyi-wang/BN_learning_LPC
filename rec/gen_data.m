% gen_data  Generate the synthetic datasets according to ODEs.
%
% 
% Is implemented by the methods used in Mendes et al. (2003) &  Soranzo et al. (2007) 
% gen_profiles.m rank_transf.m rank_transf.m netw_kinet.m are from Soranzo's software package
% mk_rnd_dag.m sample_discrete.m normalise.m  topological_sort.m were from BNT toolbox (http://www.cs.ubc.ca/~murphyk/Software/BNT/bnt.html)
% 
% Soranzo, N., G. Bianconi, and C. Altafini. 2007. 
% Comparing association network algorithms for reverse engineering of large-scale gene regulatory networks: syn-thetic versus real data. Bioinformatics 23: 1640-1647.
%
% Mendes, P., W. Sha, and K. Ye. 2003. 
% Artificial gene networks for objective compari-son of analysis algorithms. Bioinformatics 19 Suppl 2: ii122-129.
%
 

n_rec=10;
n_genes = 100;
edge_num=n_genes*5;
n_time_points = 1; % number of time points for each experiment, if 1 we have steady state experiments

knockout = true; % whether to generate data using knockout
noise_stddev = 0.01; % standard deviation of the normal distribution of the noise added to the data matrix
max_alpha = 0.005; % maximum alpha for the calculation of pAUC(ROC)


n_measurements = [10 20 50 100 200 500 1000];
output_directory = ['../data/data_' num2str(n_genes) 'genes'];
filenames_start = [num2str(n_genes) 'g_' num2str(n_rec) 'rec_' ];


if output_directory & output_directory(end) ~= '/'
	output_directory = [output_directory '/'];
end
if n_genes <= 500
	ggm_protect = 0.1;
else
	ggm_protect = 0.001;
end
n_n_measurements = numel(n_measurements);

for i = 1:n_rec
	i
    edges=0;
    while edges<edge_num
       [A, ord,edges]=mk_rnd_dag(n_genes,10,edge_num);               
    end
    
	
    dlmwrite([output_directory filenames_start 'A_' num2str(i) '.txt'], A,'\t');
	

	A = sign(A .* (rand(n_genes) - 0.5));
    V=1*ones(n_genes,1);
    theta=100*ones(n_genes);
    h= 1*ones(n_genes);

	lambda = 0.1 * rand(n_genes, 1);


	for j = 1:n_n_measurements
		j


		x = gen_profiles(A, n_measurements(j), V, theta, h, lambda, n_time_points, knockout);
		% Add Gaussian noise to x
		x = x + noise_stddev * randn(size(x));
		x(x < 0) = 0;
        dlmwrite([output_directory filenames_start 'x_' num2str(i) '_' num2str(n_measurements(j)) 'measurements.txt'], x,'\t');

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


	end
end

