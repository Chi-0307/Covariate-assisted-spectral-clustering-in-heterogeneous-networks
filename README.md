======= Descriptions ========


||SIMULATIONS||
1. Helper Functions : Contains main four functions - HetCov, HomCov, Het and Hom clustering. Other auxilliary functions are used inside them.
2. varying_r3_and_d : Contains codes to create misclustering error rate plots for Hetcov with other baselines under different setups. Varying r3 and d_cov individually and jointly 
affects the misclustering error rate of both the types of the nodes. 
3. plots: Contains plots under different setups.
4. Params: No of Types is hardcoded to T=2. No of clusters(K) = 3. Cov-dimension(R) = 3. Bernoulli variable is chosen as the nodal covariates. No. of type-1 nodes = No. of type-2 nodes, under all setup. Other parameters are specified in the code as well as paper.

||UNGA_data_analysis||
1. Helpers: codes for unga data preparation - homogeneous and heterogeneous unga networks for 1990-2000 along with covariates. It contains codes for hom,het,homcov and hetcov clusterings for each year. 
2. adj_prep: preparation of homogeneous adj.
3. unga_metrics: metrics to compute aff and lcp. used in LCP_and_AFF_plots
4. UNGA_Plots: generates cluster and world plots for 1990-2000. the arguments are cluster_map and year, got from helper functions.
5. LCP_and_AFF_plots: lcp and affiliation plots for organizations(SCO,BRICS etc) for 1990-2000.
6. ARI_Comparison: ARIs between consecutive years for each benchmark methods.
7. Boxplots_Voeten: boxplots showing distribution of ideal point estimates (Voeten) within the clusters generated from proposed method over 1990-2000.
