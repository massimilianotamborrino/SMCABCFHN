# guidedABCFHN
R package for the ABC inference for the FHN model proposed by A. Samson, M. Tamborrino and I. Tubikanec 

[1] A. Samson, M. Tamborrino, I. Tubikanec

Among all methods, we rely on the guided sequential ABC approach proposed by U.Picchini and M.Tamborrino in the paper

[2] U. Picchini, M. Tamborrino (2022), Guided sequential ABC schemes for intractable Bayesian models, arXiv, 2206.12235, https://arxiv.org/abs/2206.12235

The R-code was written by Massimiliano Tamborrino (firstname dot secondname at warwick.ac.uk) in R and Rcpp and uploaded on GitHub. The Matlab version for the guidedABC paper is available at https://github.com/umbertopicchini/guidedABC

# What can you find in the package
Here we provide the code for the classic ('standard') SMC-ABC, its optimased version ('olcml') and our guided sequential ABC approaches, in particular blocked, blockedopt, full_cond, full_cond_opt for the FHN model. Note that blocked and blocked_opt are guided SIS-ABC schemes, while full_cond and full_cond_opt are guided SMC-ABC schemes. We also accommodate the code for having both classic and "model-based" summary statistics.

The main routines are "abc_sequential_fixed_threshold" and "abc_sequential_varying_threshold", which perform guided SIS-ABC, copula-based guided SIS-ABC, guided SMC-ABC, and non-guided SMC-ABC  with either user-defined threshold vector or automatically decreasing thresholds, respectively. 

# How to install the package
* Tools/Install packages/ select the source folder
*To update The simplest way is to do it via devtools, using devtools::install_github("massimilianotamborrino/guidedABCFHN")


#### (non-copula based) Guided SIS-ABC
If you are interested in the guided SIS-ABC schemes but not based on copulas, then select "blocked", "blockedopt" or "hybrid" as samplers, and type_copula=0. The latter will allow you to sample from the Gaussian proposal function without relying on copulas. The input objecs 'nu',  'type_marginals' and 'nu_marginals' are NOT relevant for these schemes. As their entries are needed (even if not used), select them to any integer number, e.g. 5, to run the routine.

#### Copula-based guided SIS-ABC
If you are interested in the copula-based guided SIS-ABC schemes, i.e. cop-blocked, cop-blockedopt and cop-hybrid, as called in the paper, then choose "blocked", "blockedopt" or "hybrid" as samplers, and then select 
- type_copula = 1 for Gaussian copula
- type_copula = 2 for t Copula. If doing so, you would also need need to specify the degree of freedom nu>3.

You then also need to specify the desired marginal distributions, choosing from 
- type_marginals = 1 for triangular 
- type_marginals= 2 for location-scale Student't. If  doing so, you would also need to specify the degree of freedom nu_marg. 
- type_marginals= 3 for logistic
- type_marginals= 4 lognormal
- type_marginals= 5 for uniform
- type_marginals= 6 for Gaussian.

Note that input values for nu and nu_marginals are always required, even if they are only relevant for type_copula=2 and type_marginals=2, respectively. If you do not work with t-copulas or t-marginals, then set them to any number, e.g. nu=nu_marginals=5. 

#### Guided SMC-ABC
Choose the samplers "full_cond" or "full_cond_opt" for the guided SMC-ABC without or with the covariance matrix optimization, respectively. The input objects 'type_copula', 'type_marginals', 'nu' and 'nu_marginals' are NOT relevant for these schemes. As their entries are needed (even if not used), select them to any integer number to run the routine. You may choose "type_copula=0" to remind yourself that the underlying proposal function is Gaussian. 

#### Non-guided SMC-ABC
If you are interested on the non-guided SMC-ABC, please select "standard" or "olcm" for the optimized covariance matrix. The input objects 'type_copula', 'type_marginals', 'nu' and 'nu_marginals' are NOT relevant for these schemes. As their entries are needed (even if not used), select them to any integer number to run the routine. You may choose "type_copula=0" to remind yourself that the underlying proposal function is Gaussian. 


# Output files of "abc_sequential_fixed_threshold" and "abc_sequential_varying_threshold"
The output files will be saved in the specified folder, and will contain information on iteration-stage (t), attempt, type of copula, type of marginals, nu, nu_marginal, with "type_cop, type_marginals, nu, nu_marg" relevant only for the copula-based methods (so they do not have a meaning in the output files if running non-copula methods). Remark: As the output keeps track of the type of copulas and marginals, you may run all types of copulas-based samplers and marginals  in the same folder. 


Output files: The input "folder" is the name of the folder where you want your results to be stored. The whole path to access to the folder can be given as input. That folder will contain many files upon completion of a run. Here is a description:

- posterior draws: a typical file name ABCdraws_stageX_attemptY_copulaZ_margA_nuB_numargC_partN This is a matrix of d x N values containing N parameters draws (particles) for each of the d parameters to infer. The X in 'stageX' is the iteration number the draws have been sampled at. The Y in 'attemptY' is the ID of the inference run (since it is possible to run several inference runs on the same dataset, one after the other). Z, A, B and C in 'copulaZ_margA_nuB_numargC' refer to the type of copula (Z=0 means no copula method is used. See the comments above for info about Z, A, B and C). N in 'partN' denotes the desired number of accepted parameter draws, also named "particles".
- ess_stageX_attemptY_copulaZ_margA_nuB_numargC_partN: Value of the effective sample size (ESS) at the given iteration, at the given attempt etc.
- evaltime_stageX_attemptY_copulaZ_margA_nuB_numargC_partN: Number of seconds required to accept N draws at the given iteration, at the given attempt etc.
- numproposals_stageX_attemptY_copulaZ_margA_nuB_numargC_partN: Number of particles proposed (ie both the accepted and rejected), at the given iteration, attempt etc.
- weights_stageX_attemptY_copulaZ_margA_nuB_numargC_partN: Vector of N normalised importance weights associated to each of the N accepted particles at the given iteration, at the given attempt etc.
- ABCthreshold_stageX_attemptY_copulaZ_margA_nuB_numargC_partN: Value of the ABC threshold used at the given iteration, attempt etc.

