These files implement the cluster identification procedure described in
"Cluster identification using projections", by D. Peña and F.J. Prieto,
JASA, 96 (2001)

The files required by the procedure based on the projection onto
kurtosis maximization and minimization directions are:

clus_kur                       the main file (see help clus_kur)
clus_bas                       identification of significant gaps and labelling of observations
clus_grp                       a routine to analyze the convenience of merging identified clusters
ord_clus                       a routine to relabel the clusters (for housecleaning purposes)
kur_nwa                        a file that computes the orthogonal projection directions
max_kur                        a file to compute one maximization direction
min_kur                        a file to compute one minimization direction
val_kur                        a file that evaluates the projected kurtosis coefficient

The last seven files should not be used independently.

The following files implement the Jones and Sibson, Journal of the Royal Statistical Society
B (1987) criterion for the detection of clusters.

clus_js                        the main file (see help clus_js)
mcmix_js                       computation of projection directions and projected observations
max_kur                        computation of one direction maximizing the Jones and Sibson criterion
clus_grp                       as above
ord_clus                       as above

The following files have been used to conduct simulation experiments, as reported in the
paper.

sim_clus                       Conducts simulation experiments using the kurtosis projection
                               procedure (clus_kur)
sim_cljs                       Conducts simulation experiments using the Jones and Sibson
                               procedure (clus_js)
gen_clua                       Generates clusters of normal observations with outliers
gen_club                       Generates clusters from normal/uniform/Student-t observations and
                               not very ill-conditioned covariance matrices
gen_clus                       Generates clusters from normal/uniform/Student-t observations and
                               ill-conditioned covariance matrices
scp_cls                        A script file that generates the results reported in the paper
normaliz                       A file to standardize observations

Instructions to use the different files are included in the help for each file.
In particular, clus_kur expects as input a matrix having each observation as one of
its rows. The size of the matrix will be the number of observations times the
number of variables.

All files have been written using Matlab 4.2.

29/6/01