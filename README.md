# SAS-BCAR: Spatial Adaptive Selection using Binary Conditional Autoregressive Model with Application to Brain-Computer Interface



![radius_simulation](https://github.com/zikaiLin/sas-bcar/blob/main/radius_simulation.png)



## Abstarct

In medical imaging studies, statistical inferences on scalar-on-image regression are challenging due to the limited sample size and the high-dimensionality of datasets. Also, the imaging predictors often exhibit spatially heterogeneous activation patterns and the complex nonlinear associations with the response variable. To address these issues, we propose a novel prior model for Bayesian scalar-on-image regression, the Spatial Adaptive Selection using Binary Conditional Autoregressive Model (SAS-BCAR) prior.

The proposed SAS-BCAR prior employs a binary conditional autoregressive model to address spatial dependencies among feature selection indicators, effectively identifying spatially structured sparsity patterns in image datasets. Moreover, SAS-BCAR allows for adaptive feature selection to the varying spatial dependencies across different image regions, leading to a more precise and robust feature selection process for image analysis.
We have developed an efficient posterior computation algorithm for SAS-BCAR and  demonstrated the advantages of SAS-BCAR in image classification tasks with the benchmark computer vision datasets and in analysis of electroencephalography data in Brain computer interface applications.



## Code

`BasisExpansion_GP.cpp`: Core SAS-BCAR algorithm implemented in RCpp. 

`create_partition.R`: An example of creating an image prior partition,  $\Pi_k$,  using voxel variance across samples is provided. Note that this is just one example of creating an image partition, and there are other possible methods.

`ModelSelection.R`: Model selection (chain selection) using likelihood.

`SASBCAR_Fashion_MNIST.R` : An example code for simulation studies with an application to Fashion MNIST (for one replicate) is provided. Simply run the code, and the results will be saved locally. 

`update_pj_mat.cpp`: Utility function for the Gibbs sampler to update the prior partition probability map, as part of the SAS-BCAR updating scheme.



**Data analysis**



`Data Analysis/model_cluster_script_179.R`: Real data analysis script, single subject. 

`Data Analysis/character_accuracy.R`: Code for BCI experiments predictive accuracy.

`Data Analysis/plot_non_linear_GP_BCI.R`: Code for plotting the non-linear relationship detected by SAS-BCAR.