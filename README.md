# HOIs_multitrophic
Code of the article entitled "Multitrophic higher-order interactions modulate species persistence" by Lisa Buche, Ignasi Bartomeus and Oscar Godoy (https://www.biorxiv.org/content/10.1101/2021.11.18.469079v1.abstract). 

Buche, L. Bartomeus,I and Godoy, O (2023), Buchel9844/HOIs_multitrophic: Initial release (v0.1). Zenodo. 

Contact details: Lisa Buche (buchel9844@gmail.com or lbuche@student.unimelb.edu.au); Ignasi Bartomeus (nacho.bartomeus@gmail.com); Oscar Godoy (oscar.godoy@uca.edu.es)

Statement of authorship: IB and OG designed the study and collected data. LB performed modeling work and analyzed the data with inputs from OG and IG. LB wrote the first draft of the manuscript, and all authors contributed substantially to revisions.

Authorship of code: 
"code/HOI_2_fit_fecunfity_model.R" was modified from Mayfield, Margaret; Stouffer, Daniel (2017), Data from: Higher-order interactions capture unexplained complexity in diverse communities, Dryad, Dataset, https://doi.org/10.5061/dryad.3562g; original code name "ModelFit.R".
"code/HOI_3_Function_fit_alpha_r.R" was modified from Saavedra, Serguei et al. (2017), Data from: A structural approach for understanding multispecies coexistence, Dryad, Dataset, https://doi.org/10.5061/dryad.v9f5s;  original code name "toolbox_coexistence.R".
"code/HOI_5_Graphs_WPettry.R" was modified from Petry, W. and V. Lepori. (2022), wpetry/StructuralCoexistence: Initial release (v0.1). Zenodo. https://doi.org/10.5281/zenodo.7114127. 
"code/HOI_4_pesistence.R" was modified from  Saavedra, S., Medeiros, L. P., & AlAdwani, M. (2020). Structural forecasting of species persistence under changing environments. Ecology Letters. https://doi.org/10.1111/ele.13582- code on GitHub: https://github.com/MITEcology/ELE_Saavedra_etal_2020. Script used: " code/functions/lotka_volterra.R","code/functions/lv_pruning.R","code/functions/simplex_sampling.R";"code/fig2_plots.R";"code/fig2_sim.R".
Other scripts were written by LB. 


Abstract: Interactions between species are affected by the density of a third species in a multispecies context. How these higher-order interactions (HOIs) affect species persistence remains poorly understood. To explore the effect of HOIs stemming from multiple trophic layers on a plant community composition, we experimentally built a mesocosm with three plants and three pollinator species arranged in a fully nested or modified network structure. We estimated pairwise interactions among plants and between plants and pollinators, as well as HOIs initiated by a plant or a pollinator affecting plant species pairs. Using a structuralist approach, we evaluated the effects of the statistically supported HOIs on the persistence probability of each of the three competing plant species and their combinations. HOIs substantially affect the strength and sign of pairwise interactions between plant species, promoting the opportunities for multispecies communities to persist compared to a non-HOIs scenario. Eliminating experimentally a plant-pollinator interaction promotes a single-species community by modifying the per capita interaction strengths of both pairwise interactions and their HOIs. Our study provides empirical evidence of the joint importance of HOIs and network structure in determining the persistence probability of species within diverse communities.


The HOI_0_wrapper.R contains all the necessary code to compute the coefficients and the AIC for each focal species and situation. The scritps are organised according to the following: 
"code/HOI_0_wrapper.R" - calls all the following scripts in order (ordered by their name). 
"code/HOI_1_Fusion_Dataframe.R" - calls the 3 datasets and creates a list (fecundity.data) with two data frames for each focal species, with and without the link; 6 in total.
"code/HOI_2_fit_fecunfity_model.R" - contains the function running the population model.
"code/HOI_3_Function_fit_alpha_r.R" - contains the functions which extract the Omega, centroid, theta, fesibileity domaine. 
"code/HOI_4_pesistence.R" -  to make Fig.3
"code/HOI_5_Graphs_WPettry.R" - creates the diagram of Fig.3
"code/HOI_6_Procrustes.R" - creates Fig.S4

The data folder contains the following .csv (i) "pollinators.csv", (ii) "2016_no_pollinators_for_2017.csv", and (iii) "used_data2017.csv". 
"used_data2017.csv" and "2016_no_pollinators_for_2017.csv" contain the following variables:
number - link: which link of the network is present - treatment:treatment of the cage, determined by the number and identity of plants and pollinators - pot - focal: focal plant species - background_T: number of tomato individual around the focal species - background_R: number of radish individual around the focal species - background_H:  number of field bean individual around the focal species - fruit_number: number of fruit on the focal individual - exclusion - Notes - mean_weigthxseed: weight of seed for the focal individual - mean_seedsxfruit : number of seed per fruit for the focal individual - total_seeds - total_seed_weigth

"pollinators.csv" contains the pollinators sampling for the 2 years, variables are: 
year- species - treatment: treatment of the cage, determined by the number and identity of plants and pollinators -measure:the experimental parameters evaluated -value

(self-explantory variables are not explained)

Session info to run the code: R version 4.2.3 (2023-03-15)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Ventura 13.3.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] glmmADMB_0.8.3.3     cowplot_1.1.1        ggpubr_0.6.0         MuMIn_1.47.5         MASS_7.3-60          ggtern_3.4.2         ggsci_3.0.0         
 [8] RColorBrewer_1.1-3   corrplot_0.92        ggplotify_0.1.2      boot_1.3-28.1        JuliaCall_0.17.5     diffeqr_1.1.3        deSolve_1.36        
[15] scales_1.2.1         viridis_0.6.4        viridisLite_0.4.1    igraph_1.5.1         mvtnorm_1.2-2        scatterplot3d_0.3-44 vegan_2.6-4         
[22] lattice_0.20-45      permute_0.9-7        RcmdrMisc_2.7-2      sandwich_3.0-2       car_3.1-1            carData_3.0-5        ggthemes_4.2.4      
[29] plyr_1.8.8           lubridate_1.9.2      forcats_1.0.0        stringr_1.5.0        dplyr_1.1.2          purrr_1.0.1          readr_2.1.4         
[36] tidyr_1.3.0          tibble_3.2.1         ggplot2_3.4.3        tidyverse_2.0.0     

loaded via a namespace (and not attached):
 [1] colorspace_2.1-0     ggsignif_0.6.4       class_7.3-21         htmlTable_2.4.1      base64enc_0.1-3      rstudioapi_0.14      proxy_0.4-27        
 [8] hexbin_1.28.3        rstan_2.21.8         fansi_1.0.4          codetools_0.2-19     splines_4.2.3        robustbase_0.99-0    cachem_1.0.7        
[15] knitr_1.42           Formula_1.2-5        broom_1.0.4          cluster_2.1.4        latex2exp_0.9.6      compiler_4.2.3       backports_1.4.1     
[22] Matrix_1.5-3         fastmap_1.1.1        cli_3.6.0            htmltools_0.5.4      prettyunits_1.1.1    tools_4.2.3          coda_0.19-4         
[29] gtable_0.3.3         glue_1.6.2           Rcpp_1.0.10          cellranger_1.1.0     vctrs_0.6.0          nlme_3.1-162         tensorA_0.36.2      
[36] xfun_0.37            proto_1.0.0          ps_1.7.3             timechange_0.2.0     lifecycle_1.0.3      rstatix_0.7.2        DEoptimR_1.1-1      
[43] zoo_1.8-11           hms_1.1.3            parallel_4.2.3       inline_0.3.19        memoise_2.0.1        gridExtra_2.3        loo_2.5.1           
[50] StanHeaders_2.21.0-7 yulab.utils_0.0.7    rpart_4.1.19         stringi_1.7.12       nortest_1.0-4        e1071_1.7-13         checkmate_2.1.0     
[57] pkgbuild_1.4.0       R2admb_0.7.16.3      compositions_2.0-6   rlang_1.1.0          pkgconfig_2.0.3      matrixStats_0.63.0   evaluate_0.20       
[64] htmlwidgets_1.6.2    processx_3.8.0       tidyselect_1.2.0     magrittr_2.0.3       R6_2.5.1             generics_0.1.3       Hmisc_5.1-0         
[71] DBI_1.1.3            pillar_1.9.0         haven_2.5.2          foreign_0.8-84       withr_2.5.0          mgcv_1.8-42          abind_1.4-5         
[78] nnet_7.3-18          bayesm_3.1-5         crayon_1.5.2         utf8_1.2.3           tzdb_0.3.0           rmarkdown_2.20       grid_4.2.3          
[85] readxl_1.4.2         data.table_1.14.8    callr_3.7.3          digest_0.6.31        gridGraphics_0.5-1   RcppParallel_5.1.7   stats4_4.2.3        
[92] munsell_0.5.0     
