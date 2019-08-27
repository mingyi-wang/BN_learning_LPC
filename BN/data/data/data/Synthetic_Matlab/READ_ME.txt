
The synthetic Gaussian Raf-Mek-Erk pathway data are provided 
in Matlab (by Mathworks) format, but should also be accessible 
using Ocatve. See http://www.gnu.org/software/octave/ 
for further details.


The originally sampled data sets and the true gold standard network 
analysed in  Werhli et al. (2006) are also available from
http://www.bioss.ac.uk/staff/adriano/comparison/comparison.html



DATA.mat provides the five 11-by-100 data matrices from which the 
mixtures were sampled.

For i=1,...,5 

DATA{i} is the i-th data set employed by Werhli et al. (2006)

The rows correspond to the following genes/proteins:

1 raf, 2 mek, 3 plcg, 4 pip2, 5 pip3, 6 erk, 7 akt
8 pka, 9 pkc, 10 p38, 11 jnk


Each column is an iid realisation of the domain.
 
The true gold standard network from Sachs et al. (2005) is given by the 
11-by-11 matrix RAF_MEK_ERK_PATHWAY.mat

RAF_MEK_ERK_PATHWAY(i,j) = 1 
means that there is an edge pointing from the i-th node
to the j-th node, while RAF_MEK_ERK_PATHWAY(i,j) = 0 otherwise.


The five protocol Matlab structures (PROTOCOL_1.mat,...,PROTOCOL_5.mat) 
can be used to create/extract the mixture data sets from the original 
data sets in DATA.mat.


E.g.
For i=1,..,4 (corresponds to sample size m) and j=1,..,5 (corresponds to five
independent replications)

PROTOCOL_3{i}{j}

provides the indicis of the three randomly 
sampled data sets from DATA.mat, namely 
PROTOCOL_3{i}{j}.ind_1
PROTOCOL_3{i}{j}.ind_2
PROTOCOL_3{i}{j}.ind_3

and the indicis of the sampled realisations
are given by:
PROTOCOL_3{i}{j}.ind_1 
PROTOCOL_3{i}{j}.ind_2
PROTOCOL_3{i}{j}.ind_3



It should be noted that the first index i corresponds to

PROTOCOL_1.mat
i=1 -> m=30
i=2 -> m=60

PROCTOCOL_2.mat and PROTOCOL_3.mat
i=1 -> m=30
i=2 -> m=60
i=3 -> m=120
i=4 -> m=180

PROCTOCOL_4.mat and PROTOCOL_5.mat
i=1 -> m=60
i=2 -> m=120
i=3 -> m=180


The larger data sets with sample size m=480 which have been directly sampled from 
the data made freely available by Sachs et al. (2005) are given by
five Matlab structures:

DATA_480_1.mat
DATA_480_2.mat
DATA_480_3.mat
DATA_480_4.mat
DATA_480_5.mat


where each structure consists of five data sets.

E.g.: 

For i=1,...,5
DATA_480_3.mat{i} is the i-th 11-by-480 data matrix
with k=3 mixture components

The true allocation for the m=480 realisations (rows) 
is as follows:

j = 1,...,(m/k) 		-> k=1
j = (m/k)+1,...,2*(m/k) 	-> k=2
j = 2*(m/k),...,3*(m/k) = m	-> k=3



end




