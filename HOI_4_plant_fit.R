#This script calculates plant lambda and alpha matrices, and performs 
#coexistence analysis.


################################################################################
# Compute the alpha matrice and r
#################################################################################
# for each model there is two treatement, "with link" and "no link"
model.list # list of all model, attributed in the HOI_wrapper.R
model_Coefficient <- column_to_rownames(model_Coefficient,'coeff') 
#rownames(model_Coefficient)
#str(model_Coefficient)


Coefficients_alpha_r <- list() # list where the alpha and r value of each model will be stuck
# loop to extract the alpha and r value 
for(model in model.list){
for(network in c("with_link", "no_link")){
    # create data_frame for alpha and r
    matrix_alpha <- matrix(ncol=3, nrow=3)
    colnames(matrix_alpha) <- plant.id
    rownames(matrix_alpha) <- plant.id
    matrix_sd <- matrix_alpha
    matrix_direct_int <- matrix_alpha
    matrix_direct_int_sd <- matrix_alpha
    df_r <- matrix(ncol=3, nrow=1)
    rownames(df_r) <- c("r.values")
    colnames(df_r ) <- plant.id
    # create the vectors with the different HOIs
    HOI_plants_poll <- as.vector(outer(plant.id,pollinator.id, paste, sep=":"))
    imp.interaction <- c("Tomato:Lucilia_sericata","Tomato:Osmia_bicornis","Vicia:Lucilia_sericata")
    HOI_plants_poll  <- HOI_plants_poll[!HOI_plants_poll%in% imp.interaction]
    imp.interaction.no.link <- c("Raphanus:Bombus_terrestris")
    if (network=="no_link"){
        HOI_plants_poll  <- HOI_plants_poll[!HOI_plants_poll%in% imp.interaction.no.link]
    }
    
    HOI_plants <- combn(plant.id,2)
    HOI_plants <- apply(HOI_plants,2,paste,collapse=":")
    HOI_poll <- combn(pollinator.id,2)
    HOI_poll  <- apply(HOI_poll ,2,paste,collapse=":")
for(plant.focal in plant.id){
    for(plant.comp in plant.id){
matrix_alpha[plant.focal,plant.comp] <- as.numeric(model_Coefficient[plant.comp,
                                                                      paste(model,
                                                                           plant.focal,network,sep="_")]) + 
    sum(as.numeric(model_Coefficient[HOI_plants_poll[grepl(plant.comp,HOI_plants_poll)],
                          paste(model,plant.focal,network,sep="_")]),na.rm = T) +
    sum(as.numeric(model_Coefficient[HOI_plants[grepl(plant.comp,HOI_plants)],
                                     paste(model,plant.focal,network,sep="_")]),na.rm = T)
# extraction of purely pairwise interaction 
if (model == "Model.Null" | model == "Model.Mutualistic.effect" ){
    matrix_direct_int[plant.focal,plant.comp]  <- NA
} 
else{
matrix_direct_int[plant.focal,plant.comp] <- as.numeric(model_Coefficient[plant.comp,
                                                                     paste(model,
                                                                           plant.focal,network,sep="_")])
}

# extraction of standar error 
vec <- c(plant.comp,HOI_plants_poll[grepl(plant.comp,HOI_plants_poll)],
         HOI_plants[grepl(plant.comp,HOI_plants)])
mat.cov <- vcov(eval(as.name(paste("models_of",plant.focal, sep="_")))[[paste(model,network,sep="_")]])


length.vec <- c(1:length(vec))
cov <- c()
if (model == "Model.Null" | model == "Model.Mutualistic.effect" ){
    matrix_sd <- matrix(ncol=3, nrow=3)
    matrix_direct_int_sd <- matrix(ncol=3, nrow=3)
} 
else{
for(i in length.vec){
    if (!vec[i] %in% names(as.data.frame(mat.cov))) next
    cov.i <- mat.cov[vec[i],vec[i]]
    cov <- c(cov,cov.i)
    for(n in length.vec[!length.vec %in% 1:i]){
        if (!vec[n] %in% names(as.data.frame(mat.cov))) next
        cov.n <- 2*mat.cov[vec[i],vec[n]]   
        cov <- c(cov,cov.n)
    }
}
matrix_sd[plant.focal,plant.comp]  <- sqrt(sum(cov))
# extraction of the standard error of the purely pairwise interaction 
matrix_direct_int_sd[plant.focal,plant.comp] <-  sqrt(mat.cov[plant.focal,plant.comp])

}

    }
    
df_r["r.values",plant.focal] <- as.numeric(model_Coefficient["(Intercept)",
                                                             paste(model,
                                                                   plant.focal,network,sep="_")]) +
    as.numeric(model_Coefficient["year2017",
                                 paste(model,
                                       plant.focal,network,sep="_")]) + 
    sum(as.numeric(model_Coefficient[pollinator.id,paste(model,plant.focal,network,sep="_")]),na.rm = T) + 
    
    sum(as.numeric(model_Coefficient[HOI_poll,paste(model,plant.focal,network,sep="_")]),na.rm = T)

                            }
# alpha must be a matrix and r a vector

Coefficients_alpha_r[[paste(model,network, sep="_")]] <- assign(paste(model,network, sep="_"), 
                                        list(alpha=matrix_alpha, r = c(df_r["r.values",]), 
                                             sd = matrix_sd, direct_int=matrix_direct_int,
                                             direct_int_sd = matrix_direct_int_sd))
      }
}

# check 8 
  check8 <- data.frame(coeff= rownames(model_Coefficient),
             value= as.numeric(model_Coefficient$Model.HOI.plant.on.plant.and.poll.on.mutualism.and.poll.on.comp_Raphanus_with_link))
  check8 <- column_to_rownames(check8, var="coeff")
  # add all coefficient determining alpha_raphanus
  sum(check8["Raphanus","value"],check8["Raphanus:Bombus_terrestris","value"],check8["Raphanus:Lucilia_sericata","value"],
      check8["Raphanus:Osmia_bicornis","value"],check8["Raphanus:Tomato","value"],check8["Raphanus:Vicia","value"], 
      na.rm=T)
 # the value should be the same than this one:
   Coefficients_alpha_r$Model.HOI.plant.on.plant.and.poll.on.mutualism.and.poll.on.comp_with_link$alpha["Raphanus","Raphanus"]

   # add all coefficient determining r_raphanus
  sum(check8["(Intercept)","value"],check8["year2017","value"], 
      check8["Bombus_terrestris","value"],check8["Lucilia_sericata","value"],
      check8["Osmia_bicornis","value"],check8["Bombus_terrestris:Lucilia_sericata","value"],
      check8["Bombus_terrestris:Osmia_bicornis","value"],check8["Lucilia_sericata:Osmia_bicornis","value"],
      na.rm=T)
  # the value should be the same than this one:
  Coefficients_alpha_r$Model.HOI.plant.on.plant.and.poll.on.mutualism.and.poll.on.comp_with_link$r["Raphanus"]


####################################################################################################################
# Compute Omega, centroid, theta, fesibility domaine , test_feasibility_pairs and compute_overlap
##############################################################################################################
for(model in model.list){
    for(network in c("with_link", "no_link")){
alpha <- Coefficients_alpha_r[[paste(model,network, sep="_")]]$alpha
if (model == "Model.Null" | model == "Model.Mutualistic.effect" ) next #remove the null model and mutualistic model 
r <- Coefficients_alpha_r[[paste(model,network, sep="_")]]$r
#Calculate omega
Coefficients_alpha_r[[paste(model,network, sep="_")]]$Omega_fit <- 10^(Omega(norm_L2(alpha))) 
# OG: This is not correct value according to shinny app of OMEGA is 0.198 not -1.62, We lack doing the 10^ on Omega. 
r_c <- r_centroid(norm_L2(alpha)) 
Coefficients_alpha_r[[paste(model,network, sep="_")]]$theta_fit <- theta(r,r_c) 
Coefficients_alpha_r[[paste(model,network, sep="_")]]$feasibility <- test_feasibility(alpha,r)
Coefficients_alpha_r[[paste(model,network, sep="_")]]$feasibility_pair <- test_feasibility_pairs(alpha,r)
Coefficients_alpha_r[[paste(model,network, sep="_")]]$Overlap <- compute_overlap(alpha,100)

    }
}

