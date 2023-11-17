Code of the article entitled "Multitrophic higher-order interactions modulate species persistence" by Lisa Buche, Ignasi Bartomeus and Oscar Godoy (https://www.biorxiv.org/content/10.1101/2021.11.18.469079v1.abstract). 

Buche, L. Bartomeus,I and Godoy, O (2023), Buchel9844/HOIs_multitrophic: Initial release (v0.1). Zenodo. DOI: 10.5281/zenodo.8260136

Contact details: Lisa Buche (buchel9844@gmail.com or lbuche@student.unimelb.edu.au); Ignasi Bartomeus (nacho.bartomeus@gmail.com); Oscar Godoy (oscar.godoy@uca.edu.es)

Statement of authorship: IB and OG designed the study and collected data. LB performed modelling work and analyzed the data with inputs from OG and IG. LB wrote the first draft of the manuscript, and all authors contributed substantially to revisions.

Authorship of code: 
"code/HOI_2_fit_fecunfity_model.R" was modified from Mayfield, Margaret; Stouffer, Daniel (2017), Data from: Higher-order interactions capture unexplained complexity in diverse communities, Dryad, Dataset, https://doi.org/10.5061/dryad.3562g; original code name "ModelFit.R".
"code/HOI_3_Function_fit_alpha_r.R" was modified from Saavedra, Serguei et al. (2017), Data from: A structural approach for understanding multispecies coexistence, Dryad, Dataset, https://doi.org/10.5061/dryad.v9f5s;  original code name "toolbox_coexistence.R".
"code/HOI_5_Graphs_WPettry.R" was modified from Petry, W. and V. Lepori. (2022), wpetry/StructuralCoexistence: Initial release (v0.1). Zenodo. https://doi.org/10.5281/zenodo.7114127. 
"code/HOI_4_pesistence.R" was modified from  Saavedra, S., Medeiros, L. P., & AlAdwani, M. (2020). Structural forecasting of species persistence under changing environments. Ecology Letters. https://doi.org/10.1111/ele.13582- code on GitHub: https://github.com/MITEcology/ELE_Saavedra_etal_2020. Script used: " code/functions/lotka_volterra.R","code/functions/lv_pruning.R","code/functions/simplex_sampling.R";"code/fig2_plots.R";"code/fig2_sim.R".
LB wrote other scripts. 

Authorship of data: Data was collected by the authors and previously used in Bartomeus, I., Saavedra, S., Rohr, R. P., & Godoy, O. (2021). Experimental evidence for the importance of multi-trophic structure in regulating species persistence.


Abstract: Interactions between species are affected by the density of a third species in a multispecies context. How these higher-order interactions (HOIs) affect species persistence remains poorly understood. To explore the effect of HOIs stemming from multiple trophic layers on a plant community composition, we experimentally built a mesocosm with three plants and three pollinator species arranged in a fully nested or modified network structure. We estimated pairwise interactions among plants and between plants and pollinators, as well as HOIs initiated by a plant or a pollinator affecting plant species pairs. Using a structuralist approach, we evaluated the effects of the statistically supported HOIs on the persistence probability of each of the three competing plant species and their combinations. HOIs substantially affect the strength and sign of pairwise interactions between plant species, promoting the opportunities for multispecies communities to persist compared to a non-HOIs scenario. Eliminating a plant-pollinator interaction experimentally promotes a single-species community by modifying the per capita interaction strengths of both pairwise interactions and their HOIs. Our study provides empirical evidence of the joint importance of HOIs and network structure in determining the persistence probability of species within diverse communities.


The HOI_0_wrapper.R contains all the necessary code to compute the coefficients and the AIC for each focal species and situation. The scripts are organised according to the following: 
"code/HOI_0_wrapper.R" - calls all the following scripts in order (ordered by their name and run accordingly). 
"code/HOI_1_Fusion_Dataframe.R" - calls the 3 datasets and creates a list (fecundity.data) with two data frames for each focal species, with and without the link; 6 in total.
"code/HOI_2_fit_fecunfity_model.R" - contains the function running the population model.
"code/HOI_3_Function_fit_alpha_r.R" - contains the functions which extract the Omega, centroid, theta, and feasibility domaine. 
"code/HOI_4_pesistence.R" -  to make Fig.3
"code/HOI_5_Graphs_WPettry.R" - creates the diagram of Fig.3
"code/HOI_6_Procrustes.R" - creates Fig.S4

The data folder contains the following .csv (i) "pollinators.csv", (ii) "2016_no_pollinators_for_2017.csv", and (iii) "used_data2017.csv". 
"used_data2017.csv" and "2016_no_pollinators_for_2017.csv" contain the 906 focal individuals for which we counted seeds and variables:
number - link: which link of the network is present - treatment: treatment of the cage, determined by the number and identity of plants and pollinators - pot: pot number in the cage - focal: focal plant species - background_T: number of tomato individual around the focal species - background_R: number of radish individual around the focal species - background_H:  number of field bean individual around the focal species - fruit_number: number of fruit on the focal individual - Notes: notes from the observer(s) - mean_weigthxseed: weight of seed for the focal individual - mean_seedsxfruit : number of seed per fruit for the focal individual - total_seeds: total number of seeds counted for the individual - total_seed_weigth: weight of all the seeds produced by the individual

"pollinators.csv" contains the pollinators sampling for the 2 years; variables are: 
year- species - treatment: treatment of the cage, determined by the number and identity of plants and pollinators -measure: the experimental parameters evaluated -value


Session info to run the code is located in the renv.lock file. This file was created using the Renv package, function init() and snapshot(). These functions create the renv folder and the renv.lock file, which allows reproducible environments for R projects (see https://rstudio.github.io/renv/articles/renv.html).
