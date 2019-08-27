

####################################################################################
#  The low-order constraint-based algorithm and some related software and data     #
####################################################################################

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

####################################################################################
# References                                                                       #
####################################################################################

[1]Faith, J.J., B. Hayete, J.T. Thaden, I. Mogno, J. Wierzbowski, G. Cottarel, S. Kasif, J.J. Collins, and T.S. Gardner. 2007. Large-scale mapping and validation of Es-cherichia coli transcriptional regulation from a compendium of expression pro-files. PLoS Biol 5: e8.
[2]Kalisch, M. and P. Bühlmann. 2007. Estimating High-Dimensional Directed Acyclic Graphs with the PC-Algorithm. The Journal of Machine Learning Research 8: 613-636.
[3]Schafer, J. and K. Strimmer. 2005. An empirical Bayes approach to inferring large-scale gene association networks. Bioinformatics 21: 754-764.
[4]Shannon P, Markiel A, Ozier O, Baliga NS, Wang JT, Ramage D, Amin N, Schwikowski B, Ideker T. 2003. Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome Research 2003 Nov; 13(11):2498-504.
[5]Soranzo, N., G. Bianconi, and C. Altafini. 2007. Comparing association network algorithms for reverse engineering of large-scale gene regulatory networks: syn-thetic versus real data. Bioinformatics 23: 1640-1647.
[6]Werhli, A.V., M. Grzegorczyk, and D. Husmeier. 2006. Comparative evaluation of reverse engineering gene regulatory networks with relevance networks, graphical gaussian models and bayesian networks. Bioinformatics 22: 2523-2531.



