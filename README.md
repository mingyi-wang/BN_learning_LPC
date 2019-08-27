

Summary

Recently, simplified graphical modeling approaches based on low-order conditional (in-)dependence calculations have received attention because of their potential to model gene regulatory networks. Such methods are able to reconstruct large-scale gene networks with a small number of experimental measurements, at minimal computational cost. However, unlike Bayesian networks, current low-order graphical models provide no means to distinguish between cause and effect in gene regulatory relationships. To address this problem, we developed a low-order constraint-based algorithm for gene regulatory network inference. The method is capable of inferring causal directions using limited-order conditional independence tests and provides a computationally-feasible way to analyze high-dimensional datasets while maintaining high reliability. To assess the performance of our algorithm, we compared it to several existing graphical models: relevance networks; graphical Gaussian models; ARACNE; Bayesian networks; and the classical constraint-based algorithm, using realistic synthetic datasets. Furthermore, we applied our algorithm to real microarray data from Escherichia coli Affymetrix arrays and validated the results by comparison to known regulatory interactions collected in RegulonDB. The algorithm was found to be both effective and efficient at reconstructing gene regulatory networks from microarray data.

Software Package (also be accessed from http://bioinfo.noble.org/manuscript-support/lpc/)
The folders and their content

1. LPC
The implementation of the LPC algorithm 
test_sim_lpc.m show how to revoke the two major functions: learn_struct_lpc and Orientation_lpc.
Author: Mingyi Wang
Email: mwang@noble.org or mingy.wang@gmail.com.

2. rec
The scripts for generating the syentheic datasets, which were placed into the folder "data".
Some source codes are taken from Soranzo et al. (2007).

3. BN
The implementation of Bayesian network learning method, which was taken from Werhli et al. 2006.

4. data
Simulated data and E. coli real microarray data used for the tests.
The E. coli data was downloaded from http://gardnerlab.bu.edu/data/PLoS_2007/.
The data sets can be dowloaded from the supplemental web page. Please extract into this folder.

5. results
The test results.


Data Sets (including the synthetic and E.coli real data)

Supplemental file 1

Supplemental file 2 (the top 50 predicted regulatory relations in E.coli)

Supplemental file 3 (the 16,884 predicted regulatory relations in E.coli)

Citation: Inferring large-scale gene regulatory networks using a low-order constraint-based algorithm Mingyi Wang, Vagner Augusto Benedito, Patrick Xuechun Zhao and Michael Udvardi*, Mol. BioSyst., 2010, DOI:10.1039/b917571g
