################################################################################
# Introduction of the Paper
################################################################################
# Lisa Buche, Ignasi Bartomeus & Oscar godoy

################################################################################
# Import package
#################################################################################
#rm(list = ls(all = TRUE))
{install.packages("tidyverse");library(tidyverse)}
{install.packages("dplyr");library(dplyr)}
{install.packages("plyr");library(plyr)}
{install.packages("ggplot2"); library(ggplot2)}
{install.packages("ggthemes");library(ggthemes)}
{install.packages("tibble");library(tibble)}
#library(cowplot)
#library(grid)
#library(ggpubr)
# library(gtable)
{install.packages("RcmdrMisc");library(RcmdrMisc)}
{install.packages("vegan");library(vegan)} # for the procrustes fonction
{install.packages("scatterplot3d"); library(scatterplot3d)}
{install.packages("mvtnorm"); library(mvtnorm)}
{install.packages("igraph"); library(igraph)}
#{install.packages("plyr"); library(plyr)}
#if(!require(tidyr)) {install.packages("tidyr"); library(tidyr)}
{install.packages("viridis"); library(viridis)}
{install.packages("scales"); library(scales)}
{install.packages("deSolve"); library(deSolve)}
{install.packages("diffeqr"); library(diffeqr)}
{install.packages("JuliaCall"); library(JuliaCall)}
{install.packages("boot"); library(boot)}
{install.packages("ggplotify"); library(ggplotify)}
{install.packages("corrplot"); library(corrplot)}
{install.packages("RColorBrewer"); library(RColorBrewer)}
{install.packages("ggsci"); library(ggsci)}
{install.packages("ggtern"); library(ggtern)}

{install.packages("MASS"); library(MASS)} # for the neg binomial
{install.packages("MuMIn"); library(MuMIn)} # for the dredge function
library(ggpubr)
library(cowplot)
# if using "negbin-mm", you must install the "glmmADMB" packages
{install.packages("glmmADMB", 
                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                        getOption("repos")),
                type="source"); 
  library(glmmADMB)} # for the neg binomial


################################################################################
# A. Import the Data
################################################################################
#setwd("~/Eco_Bayesian")
# read in the Fusion_Dataframe, this will create a list (fecundity.data) 
# with two dataframes for each focal species, with and without the link
# check if the fecundity.data is a list of 6 elements ( 6 data_frames)
#rm(list=ls()) # remove all variable from the environment before starting the computation
source('code/HOI_1_Fusion_Dataframe.R')
# the code contains check to facilitate the identification of potential mistakes

################################################################################
# B. Run the models with the dredge and estimate the parameters
################################################################################
# ---- 1. Make the dredge function ----

# Function to calculate maximum correlation coefficient between predictor variables, retrieved from each model
max.r <- function(x){
    if(class(x)[length(class(x))] == "lm"){
        corm <- summary(x, correlation=TRUE)$correlation}
    else if(class(x) =="lmerMod"){
        corm <- cov2cor(vcov(x))}
    else if(class(x) =="lmerModLmerTest"){
        corm <- cov2cor(vcov(x))}
    else if(class(x) =="glmerMod"){
        corm <- cov2cor(vcov(x))}
    else if(class(x)=="gls"){
        corm <- summary(x)$corBeta} 
    else if(class(x)=="lme"){
        corm <- summary(x)$corFixed}
    else { print("Error: Invalid model class")}
    corm <- as.matrix(corm)
    if (length(corm)==1){
        corm <- 0
        max(abs(corm))
    } else if (length(corm)==4){
        cormf <- corm[2:nrow(corm),2:ncol(corm)]
        cormf <- 0
        max(abs(cormf))
    } else {
        cormf <- corm[2:nrow(corm),2:ncol(corm)]
        diag(cormf) <- 0
        max(abs(cormf))
    }
}

#    data <- fecundity.data[[paste(species,network,sep="_")]]
#    model <- eval(as.name(paste("models_of",species, sep="_")))[paste("Model.HOI.plant.on.plant.and.poll.on.mutualism.and.poll.on.comp",network,sep="_")]             
dredge.function <- function(data,spcs,netwk){
    
    
    # lets figure out who the potential competitors are
    not.numerical.col <- c("focal","seeds","treatment", "year")
    
    competitors <- colnames(data)[!colnames(data) %in% not.numerical.col]
    
    # remove species that are never observed to co-occur
    for(sp in names(which(colSums(data[,competitors, drop=FALSE],na.rm = T)==0))){
        data[,sp] <- NULL
    }
    # view(data)
    # lets figure out the observed competitors 
    competitors <- colnames(data)[!colnames(data) %in% not.numerical.col]
    
    ##########################################################################################################
    # Description of the different coefficient based on the ecological network
    ##########################################################################################################
   
    all.alpha <- competitors[!competitors %in% pollinator.id]
    all.gama <- competitors[!competitors %in% plant.id]
    
    all.gama.typeII <- unlist(lapply(all.gama, function(x){paste0("I(",x,"/(1 + ",x,"))")}))
    
    # For HOI on mutualistic effects-- two pollinators ----
    
    # betas between heterospecific neighbors
    if(length(all.gama)>1){
        interbetas.polinator.on.mutu <- combn(all.gama,2)
        interbetas.polinator.on.mutu.2 <- apply(interbetas.polinator.on.mutu,2,paste,collapse=":")
        
        # combine all betas together into a single variable
        HOI_poll_poll <- c(interbetas.polinator.on.mutu.2)
    }else{HOI_poll_poll <-  c() } # c() need to be replaced by intrabetas.polinator.on.mutu if soft HOI included
    
    # For HOI on Competition----
    ### fit.betas.polinator.on.comp ###
    HOI_poll_plants <- as.vector(outer(all.alpha,all.gama, paste, sep=":"))
    imp.interaction <- c("Tomato:Lucilia_sericata","Tomato:Osmia_bicornis","Vicia:Lucilia_sericata")
    imp.interaction.no.link <- c("Raphanus:Bombus_terrestris")
    HOI_poll_plants <- HOI_poll_plants[!HOI_poll_plants %in% imp.interaction]
    if (network == "no_link"){
        HOI_poll_plants<- HOI_poll_plants[!HOI_poll_plants%in% imp.interaction.no.link ]
    }
    ### fit.betas.plant.on.comp ###
    
    # betas between heterospecific neighbors
    if(length(all.alpha)>1){
        interbetas.plant.on.comp <- combn(all.alpha,2)
        interbetas.plant.on.comp <- apply(interbetas.plant.on.comp,2,paste,collapse=":")
        # combine all betas together into a single variable
        HOI_plants_plants <- c(interbetas.plant.on.comp)
    }else{HOI_plants_plants <- c()}
    
    #all.betas <- c()
    #all.betas <-c(HOI_plants_plants,all.betas.polinator.on.comp,HOI_poll_poll)
    
    ##########################################################################################################
    # Writing of the model formula 
    ##########################################################################################################   
    # start with no competition model
    base.model.formula <- "seeds ~ 1"
    
    # add a fixed effect 
    if(nlevels(as.factor(data$year))>1){
        base.model.formula <- paste0(base.model.formula,"+ year")
    }
    # drop extraneous levels that could complicate the model-fitting code
    data$year <- as.factor(data$year)
    #data$year <- droplevels(data$year)
    #levels(data$year)
    
    ##########################################################################################################
    # Model with only the Pairwise interactions
    ##########################################################################################################
    # ---- statistically select for a linear or type II response from pollinator ----
    model.formula <- paste(base.model.formula,
                           paste0(all.alpha, collapse=" + "),
                           paste0(all.gama, collapse=" + "),
                           paste0(all.gama.typeII, collapse=" + "),
                           sep=" + ")
    
    # fit the negative binomial model as described in the main text
    model.dredge.pol <- glm.nb(
        as.formula(model.formula),
        data=data,
        na.action=na.fail)
    
    Models.pol  <- dredge(model.dredge.pol,
                          rank = "AIC",extra= c(max.r,BIC),
                          fixed=c(all.alpha,"year"))
    
    result.dredge.pol <- as.data.frame(full_join(as.data.frame(coef(Models.pol)),
                                                 as.data.frame(model.sel(Models.pol)))) ##Final model selection table
    result.dredge.pol <- subset( result.dredge.pol, select=-c(year))
    
    for (n.poll  in 1:length(all.gama)){
        result.dredge.pol <-result.dredge.pol[which((is.na(result.dredge.pol[,all.gama[n.poll]])  &
                                                         !is.na(result.dredge.pol[,all.gama.typeII[n.poll]])) 
                                                     | (is.na(result.dredge.pol[,all.gama.typeII[n.poll]]) &
                                                         !is.na(result.dredge.pol[,all.gama[n.poll]]))),]
    }
    
    result.dredge.pol <- result.dredge.pol[1,]
    result.dredge.pol.response <- result.dredge.pol %>%
        dplyr::select(any_of(c(all.gama.typeII,all.gama))) 
    
    result.dredge.pol.names <- names(result.dredge.pol.response)[which(!is.na(result.dredge.pol.response[1,]))]


    
    names(result.dredge.pol)[1] <- c("Intercept")
    result.dredge.pol$species <- spcs
    result.dredge.pol$network <- netwk
  
    
    ##########################################################################################################
    # Dredge on the full model with the appropriate functional response of pollinator
    ##########################################################################################################
    
    
    model.formula <- paste(base.model.formula,
                           paste0(all.alpha, collapse=" + "),
                           paste0(result.dredge.pol.names, collapse=" + "),
                           #paste0(HOI_poll_poll, collapse=" + "),
                           paste0(HOI_poll_plants, collapse=" + "),
                           paste0(HOI_plants_plants, collapse=" + "),
                           sep=" + ")
    
    
    # fit the negative binomial model as described in the main text
    model.dredge <- glm.nb(
        as.formula(model.formula),
        data=data,
        na.action=na.fail)
    
    Allmodels  <- dredge(model.dredge,
                         rank = "AIC",extra= c(max.r,BIC),
                         fixed=c(all.alpha,"year", result.dredge.pol.names))
    ###Run dredge specifying the number of predictor variables and including the max.r function
    
    #models.selection <- get.models(Allmodels, subset = max.r<=0.6)   ##Retrieve non-collinear models (max.r <=0.6)
    # here the variables are correlated so we do not execute the previous line
    
    result.dredge <- as.data.frame(full_join(as.data.frame(coef(Allmodels)),
                                             as.data.frame(model.sel(Allmodels)))) ##Final model selection table
   
    result.dredge <- subset( result.dredge, select=-c(year))
    result.dredge.initial <- result.dredge

    result.dredge.long <- subset(result.dredge, 
                                 result.dredge$AIC < (min(result.dredge$AIC) +2)
                                 )
    result.dredge.long <- bind_rows( result.dredge[1,], 
                                     result.dredge.long[which(result.dredge.long$logLik == max(result.dredge.long$logLik)),][1,])
        
    
    names(result.dredge.long)[1] <- c("Intercept")
     
    result.dredge.long$species <- spcs
    result.dredge.long$network <- netwk
    result.dredge <- result.dredge.long[2,]
    result.dredge.long <- result.dredge.long[1:5,]
  
    # simplify the columns names of the pollinator under a type II function
    
    #names.to.change <- names(result.dredge)[which(names(result.dredge) %in% all.gama.typeII)] 
    #names(result.dredge)[which(names(result.dredge) %in% all.gama.typeII)] <-  all.gama[which(all.gama.typeII %in% names.to.change)]
    
    return(list(result.dredge.initial =result.dredge.initial,
                best.model.result = result.dredge,
                best.model.result.long = result.dredge.long,
                best.model = model.dredge, 
                best.funct.rest.result=result.dredge.pol, 
                best.funct.rest.result.model = model.dredge.pol))
    
    
}

# ---- 2. Dredge application -----

dredge.df.with.HOIs <- data.frame()
dredge.df.with.HOIs.long <- data.frame()
dredge.model.with.HOIs <- list()
dredge.df.no.HOIs <- data.frame()
dredge.model.no.HOIs <- list()
plot.pol <- list()
for (species in plant.id){
    for (network in c("with_link", "no_link")){
        dredge.result.sp <- dredge.function(fecundity.data[[paste(species,network,sep="_")]],
                                            species,network) 
        plot.pol[[paste(species,network,sep="_")]] <- dredge.result.sp[["plot"]]
        dredge.df.sp.with.HOIs <- as.data.frame(dredge.result.sp[["best.model.result"]])
        dredge.df.sp.with.HOIs.long <- as.data.frame(dredge.result.sp[["best.model.result.long"]])
        
        
        dredge.model.with.HOIs[[paste(species,network,sep="_")]] <- dredge.result.sp[["best.model"]]
       
        dredge.df.sp.no.HOIs <- as.data.frame(dredge.result.sp[["best.funct.rest.result"]])

        
        dredge.model.no.HOIs[[paste(species,network,sep="_")]] <- dredge.result.sp[["best.funct.rest.result.model"]]
        
        if(length(dredge.df.with.HOIs)==0){
            dredge.df.with.HOIs <-  dredge.df.sp.with.HOIs
            dredge.df.no.HOIs <-  dredge.df.sp.no.HOIs
            dredge.df.with.HOIs.long <- dredge.df.sp.with.HOIs.long
        }else{ 
            dredge.df.with.HOIs <- full_join(dredge.df.sp.with.HOIs,dredge.df.with.HOIs)
            dredge.df.with.HOIs.long <- full_join(dredge.df.sp.with.HOIs.long,dredge.df.with.HOIs.long)
            
        dredge.df.no.HOIs <- full_join(dredge.df.no.HOIs,dredge.df.sp.no.HOIs)}
    }
}

dredge.df.no.HOIs$model <- "no_HOIs"
dredge.df.with.HOIs$model <- "with_HOIs"
dredge.df.with.HOIs.long$model <- 'with_HOIS'
dredge.df <- bind_rows(dredge.df.with.HOIs,dredge.df.no.HOIs)
View(dredge.df)

HOI_plants_poll <- c("Bombus_terrestris:Raphanus","Bombus_terrestris:Tomato","Bombus_terrestris:Vicia",
                     "Lucilia_sericata:Raphanus","Osmia_bicornis:Vicia" ,
                     "Osmia_bicornis:Raphanus" )



HOI_plants <- c("Raphanus:Tomato","Tomato:Vicia","Raphanus:Vicia")

HOI_poll <- c("Bombus_terrestris:Osmia_bicornis","Bombus_terrestris:Lucilia_sericata",
              "Lucilia_sericata:Osmia_bicornis")

pollinator.id.typeII <- unlist(lapply(pollinator.id, function(x){paste0("I(",x,"/(1 + ",x,"))")}))


col.order <-c( "model","species","network","df","AIC","BIC","logLik","max.r","delta","weight","Intercept","year2017",
               plant.id,pollinator.id,pollinator.id.typeII,HOI_plants,HOI_plants_poll)
names(dredge.df.with.HOIs.long )[!names(dredge.df.with.HOIs.long ) %in% col.order]

#names(dredge.df)[which(!names(dredge.df) %in% col.order)]
dredge.df <- dredge.df[,col.order]
#dredge.df.with.HOIs.long <- dredge.df.with.HOIs.long[,col.order]
names(dredge.df.with.HOIs.long)[!names(dredge.df.with.HOIs.long) %in% col.order]


write.table(dredge.df,
          file.path("HOIs_Lynxc/results/dredge.result.csv"), 
          na = "NA", append = F,
          col.names =TRUE)


write.table(dredge.df.with.HOIs.long,
          file.path("HOIs_Lynxc/results/dredge.df.with.HOIs.long.csv"), 
          na = "NA", append = F,
          col.names =TRUE)


#ggsave("analysis/HOI/Figures/plot.functional.response.pdf",
#       ggarrange(plotlist = plot.pol, nrow=3, ncol=2))

# ---- 3. Full models  ----

# read in the model-fitting function - function provided with manuscript
# check if the fit.fecundity is a function 
source('HOIs_Lynxc/HOI_2_fit_fecunfity_model.R')
full_model_Coefficient <- data.frame()
full.model.of.plant.spc <- list()
for( plant.spc in plant.id){
    for( network in c("with_link", "no_link")){  
        
        full.model.of.plant.spc[[paste(plant.spc,network, sep="_")]] <- fit.fecundity.model(
            fecundity.data[[paste(plant.spc, network, sep="_")]],
            type="negbin")
        
        data.to.merge <- as.data.frame(full.model.of.plant.spc[[paste(plant.spc,network, sep="_")]]$coeff)
        data.to.merge <- remove_rownames(as.data.frame(t(data.to.merge)))
        data.to.merge$species <- plant.spc
        data.to.merge$network <- network
        data.to.merge$AIC <- full.model.of.plant.spc[[paste(plant.spc,network, sep="_")]]$model$aic
        data.to.merge$logLik <- full.model.of.plant.spc[[paste(plant.spc,network, sep="_")]]$model$twologlik/2
        data.to.merge$BIC <- BIC(full.model.of.plant.spc$Raphanus_with_link$model)
        data.to.merge$max.r <- max.r(full.model.of.plant.spc$Raphanus_with_link$model)
        data.to.merge$df <- length(full.model.of.plant.spc$Raphanus_with_link$model$coefficients)
        
        if(length(full_model_Coefficient)==0){
            full_model_Coefficient <- data.to.merge
        }else{full_model_Coefficient <- full_join(full_model_Coefficient, data.to.merge)}
        
    }
}

view(full_model_Coefficient)
full_model_Coefficient$model <- "full_model"

names(dredge.df)
names(full_model_Coefficient)[!which(names(full_model_Coefficient) %in% names(dredge.df))]

all.coeff <- full_join(dredge.df,full_model_Coefficient)

#view(all.coeff)
all.coeff.newnames <- all.coeff
all.coeff.newnames <- all.coeff.newnames %>%
    mutate(species=case_when(species == "Raphanus" ~  "Radish",
                           species == "Vicia"  ~  "Field bean" ,
                           species== "Tomato" ~  "Tomato")) %>% 
    mutate(network =case_when(network== "no_link" ~  "Link physically prevented",
                              network == "with_link"  ~  "Nested structure"
    )) %>%
    mutate(model =case_when(model == "no_HOIs" ~  "Pairwise interactions",
                            model == "with_HOIs"  ~  "Selected HOIs",
                            model == "full_model" ~  "All HOIs"
    ))

all.coeff.newnames <- arrange(transform(all.coeff.newnames,
                               model=factor(model,level=c("Pairwise interactions","Selected HOIs","All HOIs"))),
                     model)

all.coeff.newnames <- arrange(transform(all.coeff.newnames,
                               network=factor(network,level=c("Nested structure","Link physically prevented"))),
                     network)

all.coeff.newnames <- arrange(transform(all.coeff.newnames,
                               species=factor(species,level=c("Radish","Field bean","Tomato"))),
                     species)
all.coeff.newnames <- all.coeff.newnames%>% mutate_if(is.numeric, round, digits=3)
write.table(all.coeff,
          file.path("HOIs_Lynxc/results",
                    "all.coeff.csv"), 
          na = "NA", append = F,
          col.names =TRUE)

#all.coeff <- read.csv("HOIs_Lynxc/results/all.coeff.csv",sep=' ')

# ---- graph of both type of response ----
data.raphanus <- fecundity.data[[paste("Raphanus", "with_link", sep="_")]]
raphanus.model.type <- fit.fecundity.model(data.raphanus,
    type="negbin")

data.raphanus$fitted.value.type.II <-  raphanus.model.type$model.type.II$fitted.values
data.raphanus$fitted.value.type.I <-  raphanus.model.type$model.type.I$fitted.values
data.raphanus$visit <- exp((data.raphanus$Bombus_terrestris-1)*(maxpol.log/maxplant)) + exp((data.raphanus$Lucilia_sericata -1)*(maxpol.log/maxplant)) + exp((data.raphanus$Osmia_bicornis-1)*(maxpol.log/maxplant))

write.table(data.raphanus,
          file.path("HOIs_Lynxc/results",
                    "data.raphanus.csv"), 
          na = "NA", append = F,
          col.names =TRUE)

plot.function.type <- ggplot(data=data.raphanus, aes(y=seeds)) +
    
    #geom_abline(intercept = 0, slope = 1) +
    #geom_abline(intercept = 0, slope = summary(model.dredge.pol)$dispersion,
    #            color="green") +
    stat_smooth(aes(x=visit), method="loess", se = FALSE,color="red") +
    xlab("Visitation rate") + ylab(paste("seeds of",spcs)) +
    stat_smooth(aes(x=fitted.value.type.I),method="loess", se = FALSE,color="green") +
    stat_smooth(aes(x=fitted.value.type.II),method="loess", se = FALSE,color="blue") +

    theme_bw() 
ggsave("HOIs_Lynxc/results/plot.function.type.pdf",
       plot.function.type)


##########################################################################################################
# C. Dredge interaction matrix
##########################################################################################################

# ---- 1. Compute the alpha matrice and r ----
Coefficients_alpha_r <- list()

for(HOIs in  c("with_HOIs","no_HOIs","full_model")){
for(network in  c("with_link", "no_link")){
    # create data_frame for alpha and r
    matrix_alpha <- matrix(ncol=3, nrow=3)
    colnames(matrix_alpha) <- plant.id
    rownames(matrix_alpha) <- plant.id
    matrix_sd <- matrix_alpha
    matrix_direct_int <- matrix_alpha
    matrix_HOIs <- matrix_alpha
    matrix_direct_int_sd <- matrix_alpha
    df_r <- matrix(ncol=3, nrow=1)
    rownames(df_r) <- c("r.values")
    colnames(df_r ) <- plant.id
    
    
    for(plant.focal in plant.id){
        for(plant.comp in plant.id){
            matrix_alpha[plant.focal,plant.comp] <- as.numeric(all.coeff[which(all.coeff$species==plant.focal &
                                                                                   all.coeff$network == network &
                                                                                   all.coeff$model == HOIs ),
                                                                         plant.comp]) + 
                sum(as.numeric(all.coeff[which(all.coeff$species==plant.focal&
                                                   all.coeff$network == network&
                                                   all.coeff$model == HOIs),
                                         HOI_plants_poll[grepl(plant.comp,HOI_plants_poll)]]),na.rm = T) +
                sum(as.numeric(all.coeff[which(all.coeff$species==plant.focal &
                                                   all.coeff$network == network &
                                                   all.coeff$model == HOIs),
                                         HOI_plants[grepl(plant.comp,HOI_plants)]]),na.rm = T)
            # direct interaction matrix
            matrix_direct_int[plant.focal,plant.comp] <- as.numeric(all.coeff[which(all.coeff$species==plant.focal &
                                                                                   all.coeff$network == network &
                                                                                   all.coeff$model == HOIs ),
                                                                         plant.comp])
            
            matrix_HOIs[plant.focal,plant.comp] <- sum(as.numeric(all.coeff[which(all.coeff$species==plant.focal&
                                                   all.coeff$network == network&
                                                   all.coeff$model == HOIs),
                                         HOI_plants_poll[grepl(plant.comp,HOI_plants_poll)]]),na.rm = T) +
                sum(as.numeric(all.coeff[which(all.coeff$species==plant.focal &
                                                   all.coeff$network == network &
                                                   all.coeff$model == HOIs),
                                         HOI_plants[grepl(plant.comp,HOI_plants)]]),na.rm = T)
            
            # extraction of standard error 
            vec <- c(plant.comp,HOI_plants_poll[grepl(plant.comp,HOI_plants_poll)],
                     HOI_plants[grepl(plant.comp,HOI_plants)])

            if (HOIs == "with_HOIs"){ mat.cov <- stats::vcov(full.model.of.plant.spc[[paste(plant.comp,network,sep="_")]]$model)}
            if (HOIs == "no_HOIs"){ mat.cov <- stats::vcov(dredge.model.no.HOIs[[paste(plant.comp,network,sep="_")]])}
            if (HOIs == "full_model"){ mat.cov <- stats::vcov( full.model.of.plant.spc[[paste(plant.spc,network, sep="_")]][[1]])}
                    
            length.vec <- c(1:length(vec))
            cov <- c()
            for(i in length.vec){
                if (!vec[i] %in% names(as.data.frame(mat.cov))) next
                cov.i <- mat.cov[vec[i],vec[i]]
                cov <- c(cov,cov.i)
                for(n in length.vec[!length.vec %in% 1:i]){
                    if (!vec[n] %in% names(as.data.frame(mat.cov))) next
                    cov.n <- 2*mat.cov[vec[i],vec[n]]   
                    cov <- c(cov,cov.n)
                }
                
                matrix_sd[plant.focal,plant.comp]  <- sqrt(sum(cov))
                # extraction of the standard error of the purely pairwise interaction 
                matrix_direct_int_sd[plant.focal,plant.comp] <-  sqrt(mat.cov[plant.focal,plant.comp])
                
            }
            
        }
        
        df_r["r.values",plant.focal] <- as.numeric(all.coeff[which(all.coeff$species==plant.focal &
                                                                       all.coeff$network == network &
                                                                       all.coeff$model == HOIs),
                                                             "Intercept"]) +
            as.numeric(all.coeff[which(all.coeff$species==plant.focal &
                                           all.coeff$network == network &
                                           all.coeff$model == HOIs),
                                 "year2017"]) + 
            sum(as.numeric(all.coeff[which(all.coeff$species==plant.focal &
                                               all.coeff$network == network &
                                               all.coeff$model == HOIs),
                                     pollinator.id]),na.rm = T) + 
            sum(as.numeric(all.coeff[which(all.coeff$species==plant.focal &
                                               all.coeff$network == network &
                                               all.coeff$model == HOIs),
                                     pollinator.id.typeII]),na.rm = T) 
        
    }
    # alpha must be a matrix and r a vector
    
    Coefficients_alpha_r[[paste(network,HOIs,sep='_')]] <- assign(paste(network,HOIs,sep='_'), 
                                                     list(alpha=matrix_alpha, r = c(df_r["r.values",]), 
                                                          sd = matrix_sd, direct_int=matrix_direct_int,
                                                          direct_int_sd = matrix_direct_int_sd,
                                                          matrix_HOIs = matrix_HOIs))
}
}

# ---- 2. Compute Omega, centroid, theta, fesibileity domaine , test_feasibility_pairs and compute_overlap ----

source('HOIs_Lynxc/HOI_3_Function_fit_alpha_r.R')
for(HOIs in  c("with_HOIs","no_HOIs","full_model")){
    for(network in  c("with_link", "no_link")){
        alpha <- Coefficients_alpha_r[[paste(network,HOIs,sep='_')]]$alpha
        r <- Coefficients_alpha_r[[paste(network,HOIs,sep='_')]]$r
        #Calculate omega
        Coefficients_alpha_r[[paste(network,HOIs,sep='_')]]$Omega_fit <- 10^(Omega(norm_L2(alpha))) 
        # OG: This is not correct value according to shinny app of OMEGA is 0.198 not -1.62, We lack doing the 10^ on Omega. 
        r_c <- r_centroid(norm_L2(alpha)) 
        Coefficients_alpha_r[[paste(network,HOIs,sep='_')]]$theta_fit <- theta(r,r_c) 
        Coefficients_alpha_r[[paste(network,HOIs,sep='_')]]$feasibility <- test_feasibility(alpha,r)
        Coefficients_alpha_r[[paste(network,HOIs,sep='_')]]$feasibility_pair <- test_feasibility_pairs(alpha,r)
        Coefficients_alpha_r[[paste(network,HOIs,sep='_')]]$Overlap <- compute_overlap(alpha,100)
        
    }
}

# ---- 3. save data in RData file ----
save(     Coefficients_alpha_r, 
          file="HOIs_Lynxc/results/Coefficients_alpha_r.RData")
#load("HOIs_Lynxc/results/Coefficients_alpha_r.RData")


################################################################################
# D. Persistence of species 
################################################################################
# the function to create the graphs has been developped by W. Petry
source('HOIs_Lynxc/HOI_4_pesistence.R')


################################################################################
# E. Create the complementary graphs
################################################################################

# ----1. Interaction network plot ----

# the function to create the graphs has been developped by W. Petry
source('HOIs_Lynxc/HOI_5_Graphs_WPettry.R')

Interaction_network <- list()
for(network in  c("with_link", "no_link")){
    for(HOIs in  c("with_HOIs","no_HOIs","full_model")){
        #g <- list()
        mk_graph_3sp(Coefficients_alpha_r[[paste(network,HOIs,sep='_')]]$alpha, 
                     Coefficients_alpha_r[[paste(network,HOIs,sep='_')]]$r)
        Interaction_network[[paste(network,HOIs,sep='_')]] <- last_plot()
        ggsave(paste0("HOIs_Lynxc/results/","network_interaction_", network,HOIs,".pdf"),
               last_plot(),
               width = 2.6,height = 2.6, units = "in")
    }
}

plot.sp.interaction.with_link <- ggarrange(plotlist=Interaction_network[1:3],ncol=3,nrow=1)

ggsave("HOIs_Lynxc/results/plot.sp.interaction.with_link.pdf",
       plot.sp.interaction.with_link, height = 3, width = 9, units = "in")

plot.sp.interaction.network <- ggarrange(plotlist=Interaction_network[4:6],ncol=3,nrow=1)
ggsave("HOIs_Lynxc/results/plot.sp.interaction.no_link.pdf",
       plot.sp.interaction.network,  height = 3, width = 9, units = "in")

# ----2. Graph of the effect of direct int vs HOIs -----

graph.dredge.interaction <- data.frame()
graph.dredge.interaction.sd <- data.frame()
for(HOIs in  c("with_HOIs","no_HOIs","full_model")){
    for(network in  c("with_link", "no_link")){
    for(plant.focal in plant.id){
        for(plant.comp in plant.id){
            graph.dredge.interaction.sp <- data.frame( network=network, model= HOIs, focal=plant.focal, 
                                                       comp=plant.comp, alpha=NA, beta.pol = NA, beta.plant=NA)
            
            graph.dredge.interaction.sp$alpha <-as.numeric(all.coeff[which(all.coeff$species==plant.focal &
                                                                               all.coeff$network == network &
                                                                               all.coeff$model == HOIs),
                                                                     plant.comp])
            
            graph.dredge.interaction.sp$beta.pol <- sum(as.numeric(all.coeff[which(all.coeff$species==plant.focal &
                                                                                       all.coeff$network == network &
                                                                                       all.coeff$model == HOIs),
                                                                             HOI_plants_poll[grepl(plant.comp,HOI_plants_poll)]]),na.rm = T)
            
            graph.dredge.interaction.sp$beta.plant <-   sum(as.numeric(all.coeff[which(all.coeff$species==plant.focal &
                                                                                           all.coeff$network == network &
                                                                                           all.coeff$model == HOIs),
                                                                                 HOI_plants[grepl(plant.comp,HOI_plants)]]),na.rm = T)
            
            graph.dredge.interaction.sp <- gather(graph.dredge.interaction.sp,alpha,beta.pol,beta.plant,
                                                  key="effect", value="value" )
            graph.dredge.interaction.sp  <- drop_na(graph.dredge.interaction.sp  ,"value")
            graph.dredge.interaction  <- bind_rows(graph.dredge.interaction,graph.dredge.interaction.sp)
            
            # extraction of standar error 
            graph.dredge.interaction.sp <- data.frame( network=network, focal=plant.focal, model=HOIs,
                                                       comp=plant.comp, alpha=NA, beta.pol = NA, beta.plant=NA)
             
            if (HOIs == "with_HOIs"){ mat.cov <- stats::vcov(dredge.model.with.HOIs[[paste(plant.comp,network,sep="_")]])}
            if (HOIs == "no_HOIs"){ mat.cov <- stats::vcov(dredge.model.no.HOIs[[paste(plant.comp,network,sep="_")]])}
            if (HOIs == "full_model"){ mat.cov <- stats::vcov(full.model.of.plant.spc[[paste(plant.spc,network, sep="_")]][[1]])}
            
            mat.names <- as.vector(names(as.data.frame(mat.cov)))
            
            vec <- c(mat.names[grepl(plant.comp,
                                         mat.names)])
            vec <-grep(plant.comp,vec,value=TRUE)
                
                graph.dredge.interaction.sp$alpha <- sqrt(as.numeric(mat.cov[plant.focal, plant.comp]))
                
                # note on Sd, var, and cov
                # SD(X) = sqrt(Var(X))
                # v(X,X)=Var(X)
                # Var(X,Y) = Var(X) + Var(Y) + 2 cov(X,Y) if X and Y are uncorrelated (i.e. Cov(X,Y)/ sqrt(Var(X) +Var(Y)) =/= 0)
                # Var(X,Y) = Var(X)
          
            
            if( length(vec) > 1 ) {
            vec.pol <- grep(paste(pollinator.id,collapse="|"),vec,value=TRUE)
            vec.plant <- grep(paste(plant.id[which(!plant.id== plant.comp)],collapse="|"),vec,value=TRUE)
            
            if (length(vec.plant) >0) {
            length.vec <- c(1:length(vec.plant))
            cov <- c()
            for(i in length.vec){
                if (!vec.plant[i] %in% names(as.data.frame(mat.cov))) next
                for(n in length.vec[!length.vec %in% 1:i]){
                    if (!vec.plant[n] %in% names(as.data.frame(mat.cov))) next
                    cov <- 2*mat.cov[vec.plant[i],vec.plant[n]]   
                }
            }
            }
            graph.dredge.interaction.sp$beta.plant <- sqrt(sum(abs(as.numeric(cov))))
           
            if (length(vec.pol) >0) {
            length.vec <- c(1:length(vec.pol))
            cov <- c()
            for(i in length.vec){
                if (!vec.plant[i] %in% names(as.data.frame(mat.cov))) next
                for(n in length.vec[!length.vec %in% 1:i]){
                    if (!vec.pol[n] %in% names(as.data.frame(mat.cov))) next
                    cov <- 2*mat.cov[vec.pol[i],vec.pol[n]]   
                }
            }
            }
            graph.dredge.interaction.sp$beta.pol <- sqrt(sum(abs(as.numeric(cov))))
            }
            
            graph.dredge.interaction.sp <- gather(graph.dredge.interaction.sp,alpha,beta.pol,beta.plant,
                                                  key="effect", value="value" )
            graph.dredge.interaction.sp  <- drop_na(graph.dredge.interaction.sp  ,"value")
            graph.dredge.interaction.sd  <- bind_rows(graph.dredge.interaction.sd,graph.dredge.interaction.sp)
        }
    }
}
}
graph.dredge.interaction <- left_join(graph.dredge.interaction,
                                      graph.dredge.interaction.sd,
                                      suffix = c("", ".sd"),
                                      keep = FALSE,
                                      by = names(graph.dredge.interaction)[which(!names(graph.dredge.interaction) %in% "value")])




graph.dredge.interaction <- graph.dredge.interaction %>%
    mutate(focal=case_when(focal == "Raphanus" ~  "Radish",
                           focal == "Vicia"  ~  "Field bean" ,
                           focal== "Tomato" ~  "Tomato")) %>% 
    mutate(comp =case_when(comp == "Raphanus" ~  "Radish",
                           comp == "Vicia"  ~  "Field bean" ,
                           comp== "Tomato" ~  "Tomato")) %>%
    mutate(network =case_when(network== "no_link" ~  "Link physically prevented",
                              network == "with_link"  ~  "Nested structure"
    )) %>%
    mutate(model =case_when(model == "no_HOIs" ~  "Pairwise interactions",
                            model == "with_HOIs"  ~  "Selected HOIs",
                            model == "full_model" ~  "All HOIs"
    ))

graph.dredge.interaction<- arrange(transform(graph.dredge.interaction,
                                             focal=factor(focal,level=c("Radish","Field bean","Tomato"))),
                                   focal)


graph.dredge.interaction <- unite(graph.dredge.interaction,focal,comp,col="interaction", sep="-",remove=F)
 
graph.dredge.interaction <- arrange(transform(graph.dredge.interaction,
                                             focal=factor(focal,level=c("Radish","Field bean","Tomato"))),
                                   focal)

graph.dredge.interaction <- arrange(transform(graph.dredge.interaction,
                                              model=factor(model,level=c("Pairwise interactions",
                                                                         "Selected HOIs","All HOIs"))),
                                    model)

graph.dredge.interaction <- drop_na(graph.dredge.interaction,value.sd)

view(graph.dredge.interaction)
write.table(graph.dredge.interaction,
          "HOIs_Lynxc/results/graph.dredge.interaction.csv")
# ---- 2.a. Test for diff between sd and HOIs effect ----
graph.dredge.interaction.alpha <- graph.dredge.interaction[graph.dredge.interaction$effect=="alpha",]
graph.dredge.interaction.beta.pol <- graph.dredge.interaction[graph.dredge.interaction$effect=="beta.pol",]
graph.dredge.interaction.beta.plant <- graph.dredge.interaction[graph.dredge.interaction$effect=="beta.plant",]

graph.dredge.interaction.effect <- left_join(graph.dredge.interaction.alpha,
                                             graph.dredge.interaction.beta.pol,
                                             by=c("network","model","interaction","focal","comp"),
                                             suffix=c(".alpha",".beta.pol")) 

graph.dredge.interaction.effect <- left_join(graph.dredge.interaction.effect ,
                                             graph.dredge.interaction.beta.plant,
                                             by=c("network","model","interaction","focal","comp"),
                                             suffix=c("",".beta.plant")) 


view(graph.dredge.interaction.effect)
graph.dredge.interaction.effect$sd.vs.beta.pol <- NA
graph.dredge.interaction.effect$sd.vs.beta.pol[which(abs(graph.dredge.interaction.effect$value.sd.alpha) <
                                                         abs(graph.dredge.interaction.effect$value.beta.pol))] <- 1
graph.dredge.interaction.effect$sd.vs.beta.pol[which(abs(graph.dredge.interaction.effect$value.sd.alpha) >
                                                         abs(graph.dredge.interaction.effect$value.beta.pol))] <- 0
graph.dredge.interaction.effect$sd.vs.beta.pol[which(graph.dredge.interaction.effect$value.beta.pol ==0)] <- NA

graph.dredge.interaction.effect$sd.vs.beta.plant <- NA
graph.dredge.interaction.effect$sd.vs.beta.plant[which(abs(graph.dredge.interaction.effect$value.sd.alpha) >
                                                         abs(graph.dredge.interaction.effect$value))] <- 0
graph.dredge.interaction.effect$sd.vs.beta.plant[which(abs(graph.dredge.interaction.effect$value.sd.alpha) <
                                                           abs(graph.dredge.interaction.effect$value))] <- 1
graph.dredge.interaction.effect$sd.vs.beta.plant[which(graph.dredge.interaction.effect$value ==0)] <- NA

graph.dredge.interaction.effect$effect.beta.pol <- NA
graph.dredge.interaction.effect$effect.beta.pol[which(graph.dredge.interaction.effect$value.beta.pol>0)] <- 1
graph.dredge.interaction.effect$effect.beta.pol[which(graph.dredge.interaction.effect$value.beta.pol<0)] <- 0
graph.dredge.interaction.effect$effect.beta.plant <- NA
graph.dredge.interaction.effect$effect.beta.plant[which(graph.dredge.interaction.effect$value>0)] <- 1
graph.dredge.interaction.effect$effect.beta.plant[which(graph.dredge.interaction.effect$value<0)] <- 0


df.effect.HOIs <- data.frame()
for(HOIs in  c("Pairwise interactions","Selected HOIs","All HOIs")){
    for(network in  c( "Link physically prevented","Nested structure")){
        df.effect.HOIs.sp <- data.frame(model= HOIs, network=network)
        
        y  <- graph.dredge.interaction.effect[which(graph.dredge.interaction.effect$network == network &
                                                    graph.dredge.interaction.effect$model == HOIs),]
        
        for(test in c("effect.beta.pol","effect.beta.plant","sd.vs.beta.plant","sd.vs.beta.pol")){
            if ( all(is.na(y[,test]))){
                df.effect.HOIs.sp[test] <- 0
            }else{ 
                
                df.effect.HOIs.sp[,test] <- round(sum(y[,test],na.rm = T)/nrow(y[!is.na(y[,test]),]),3)
            }
                
        }
        df.effect.HOIs.sp[,"Number of plant HOIs"] <- nrow(y[which(!y$value==0),])
        df.effect.HOIs.sp[,"Number of pol HOIs"] <- nrow(y[which(!y$value.beta.pol==0),])
        
        if( HOIs == "Pairwise interactions" &
            network=="Link physically prevented"){df.effect.HOIs <- df.effect.HOIs.sp 
            }else{ df.effect.HOIs <- full_join(df.effect.HOIs.sp,
                                df.effect.HOIs)
        }

    }
}
view(df.effect.HOIs)
write.table(df.effect.HOIs,
          "HOIs_Lynxc/results/df.effect.HOIs.csv")


# ---- 2.b. Final values of alpha matrix ----
coefficient.dredge_comparison <- data.frame()
coefficient.dredge_comparison.sd  <- data.frame()

for ( network in   c("with_link","no_link")){
    for(HOIs in  c("with_HOIs","no_HOIs","full_model")){
    alpha_coeff <-   as.data.frame((Coefficients_alpha_r[[paste(network,HOIs,sep="_")]]$alpha))
    alpha_coeff$network <- network
    alpha_coeff$model <- HOIs
    alpha_coeff  <- rownames_to_column(alpha_coeff ,var='focal')
    alpha_coeff.sd <-   as.data.frame((Coefficients_alpha_r[[paste(network,HOIs,sep="_")]]$sd))
    alpha_coeff.sd$network <- network
    alpha_coeff.sd$model <- HOIs
    alpha_coeff.sd  <- rownames_to_column(alpha_coeff.sd ,var='focal')
    
    coefficient.dredge_comparison <- rbind(coefficient.dredge_comparison,alpha_coeff)
    coefficient.dredge_comparison.sd <- rbind(coefficient.dredge_comparison.sd,alpha_coeff.sd)
    
    }
}




coefficient.dredge_comparison <- gather(coefficient.dredge_comparison , all_of(plant.id), 
                                        key = "competitor", value = "alpha")
coefficient.dredge_comparison.sd <- gather(coefficient.dredge_comparison.sd ,all_of(plant.id), 
                                           key = "competitor", value = "alpha.sd")

coefficient.dredge_comparison.full <- merge(coefficient.dredge_comparison,
                                            coefficient.dredge_comparison.sd)

coefficient.dredge_comparison.full <- coefficient.dredge_comparison.full %>%
    mutate(competitor.2 =case_when(competitor == "Raphanus" ~  "Radish",
                                   competitor == "Vicia"  ~  "Field bean",
                                   competitor == "Tomato" ~  "Tomato"
    )) %>%
    mutate(focal.2 =case_when(focal == "Raphanus" ~  "Radish",
                              focal == "Vicia"  ~  "Field bean",
                              focal == "Tomato" ~  "Tomato"
    ))%>%
    mutate(network =case_when(network== "no_link" ~  "Link physically prevented",
                              network == "with_link"  ~  "Nested structure"
    )) %>%
    mutate(model =case_when(model == "no_HOIs" ~  "Pairwise interactions",
                              model == "with_HOIs"  ~  "Selected HOIs",
                              model == "full_model" ~  "All HOIs"
    ))



coefficient.dredge_comparison.full <- unite(coefficient.dredge_comparison.full,focal.2, competitor.2,
                                            col="interaction", sep = "-", remove = F)
coefficient.dredge_comparison.full$sign.int <- "Net competition"
coefficient.dredge_comparison.full[which(coefficient.dredge_comparison.full$alpha > 0),"sign.int"] <- "Net facilitation"


write.table(coefficient.dredge_comparison.full,
          "HOIs_Lynxc/results/coefficient.dredge_comparison.full.csv")

# ---- 2.c. Plot ----

graph.dredge.interaction$interaction <- as.factor(graph.dredge.interaction$interaction)
grid.name <- c("Link physically prevented","Nested structure","Direct interactions","Selected HOIs","All HOIs")
names(grid.name)<- c("Link physically prevented","Nested structure","Direct interactions","Selected HOIs","All HOIs")

graph.dredge.interaction <- arrange(transform(graph.dredge.interaction,
                                              model=factor(model,level=c("Pairwise interactions","Selected HOIs","All HOIs"))),
                                    model)
coefficient.dredge_comparison.full <- arrange(transform(coefficient.dredge_comparison.full,
                                                        model=factor(model,level=c("Pairwise interactions","Selected HOIs","All HOIs"))),
                                              model)

graph.dredge.interaction <- arrange(transform(graph.dredge.interaction,
                                              network=factor(network,level=c("Nested structure","Link physically prevented"))),
                                    network)
coefficient.dredge_comparison.full <- arrange(transform(coefficient.dredge_comparison.full,
                                                        network=factor(network,level=c("Nested structure","Link physically prevented"))),
                                              network)
view(graph.dredge.interaction)
effect.dredge_interaction <- ggplot() + geom_bar(stat="identity",
                                                 alpha= 0.8,
                                                 inherit.aes = F,
                                                 data=graph.dredge.interaction,
                                                 aes(y=value, x=interaction, fill=effect)) +
    geom_point(data=coefficient.dredge_comparison.full, alpha= 0.8,
               aes(y=alpha, x=interaction, shape= sign.int)) +
    geom_linerange(data=coefficient.dredge_comparison.full,alpha= 0.8,
                   aes(y=alpha, x=interaction, ymin= alpha - alpha.sd,ymax= alpha + alpha.sd )) + 
    facet_wrap(network~model,
               labeller = labeller(.multi_line = FALSE)) +
    ylab("Interaction strength and sign") + 
    xlab("Identity of the receiver - transmitter")+
    geom_hline(yintercept=0, linetype="dashed", color = "black") + 
    scale_shape_discrete(name="") +
    scale_fill_manual(name="",
                      labels=c("Pairwise interaction","HOIs plant","HOIs pollinator"),
                      values=c("#E64B35FF","#00A087FF","#4DBBD5FF"))+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 68,vjust = 0.5),
          strip.background=element_rect(fill="white", color="white"),
          strip.text = element_text(size=10)) 
ggsave("HOIs_Lynxc/results/effect.dredge_interaction.pdf",
       effect.dredge_interaction,height = 8.27, width = 15.69, units = "in")


graph.dredge.interaction <- arrange(transform(graph.dredge.interaction,
                                      effect=factor(effect,level=c("alpha","beta.pol","beta.plant"))),
                            effect)
graph.dredge.interaction 
Overall_effect_HOIs <- ggplot(data=graph.dredge.interaction,
                              aes(y=abs(value), x=focal,fill=effect)) + 
    geom_boxplot(position = position_dodge(width = 1)) +
    geom_pointrange(aes(color=effect, ymax=(abs(value) + abs(value.sd)),
                            ymin=(abs(value) - abs(value.sd))
                      ),
                  position = position_dodge(width = 1)) +
    theme_bw()   +  
    scale_color_manual(name="",
                      values=c("#E64B35FF","#4DBBD5FF","#00A087FF"),
                      labels=c("Pairwise interaction", "HOIs pollinator","HOIs plant")) + 
    scale_fill_manual(name="",
                     values=alpha(c("#E64B35FF","#4DBBD5FF","#00A087FF"),0.5),
                     labels=c("Pairwise interaction", "HOIs pollinator","HOIs plant")) + 
  ylab("absolute value of the effect") + 
    theme(axis.text.x = element_text(face = "italic"),legend.position = "bottom")



ggsave("HOIs_Lynxc/results/Overall_effect_HOIs.pdf",
       Overall_effect_HOIs, height = 7.27, width = 8.69, units = "in")

# ----3. Graph of the persistence ----
df_prob_surv <- data.frame()
for ( network in   c("with_link","no_link")){
    for(HOIs in  c("with_HOIs","no_HOIs","full_model")){
        df_surv_i <-read.csv(paste0("HOIs_Lynxc/results/probability_surv_sp_", paste(network,HOIs,sep="_"), 
                                    ".csv", sep = ""))
        df_surv_i$model <-  HOIs 
        df_surv_i$network <-network
        df_prob_surv <- bind_rows(df_prob_surv,df_surv_i )
    }
}
df_prob_surv<- arrange(transform(df_prob_surv,
                                 network =factor(network ,level=c("with_link","no_link"))),
                       network )

df_prob_surv<- arrange(transform(df_prob_surv,
                                 model=factor(model,level=c("no_HOIs","with_HOIs","full_model"))),
                       model)

write.table(df_prob_surv,
          "HOIs_Lynxc/results/df_prob_surv.csv")

#df_prob_surv <- subset(df_prob_surv,model== model.name[1])
df_prob_surv$network <- as.factor(df_prob_surv$network)
#remove.packages("ggtern")
library("ggplot2")
df_prob_surv  <- df_prob_surv[which(df_prob_surv$model == "no_HOIs" |
                                        df_prob_surv$model == "with_HOIs"),]

df_prob_surv <- arrange(transform(df_prob_surv,
                                  network=factor( network,level=c("with_link","no_link"))),
                        network)

plot_prob_surv <- ggplot(df_prob_surv, aes(y=frac, x= sp, fill=model,
                                           color=model)) + geom_boxplot(position=position_dodge(0.4),
                                                                        alpha=0.2) +
    theme_bw() +
    facet_grid( .~network,labeller = as_labeller(c("with_link"="(a) Nested structure",
                                                   "no_link"="(b) Link physically prevented"))
    ) +
    ylab("probability of persistence") + xlab("Focal species") +
    scale_fill_manual(values=c("#FFA500","#DA70D6"),
                      labels=c("Pairwise interactions","Selected HOIs"),
                      name="") +
    scale_color_manual(values=c("#FFA500","#DA70D6")) +
    theme(strip.background=element_rect(fill="white", color="white")) + 
    guides(color= "none",fill = guide_legend(override.aes = list(alpha=0.6)))


ggsave("HOIs_Lynxc/results/plot_prob_surv.pdf",
       plot_prob_surv, height = 3.68, width = 8.51, units = "in")






# ----4. Graph of the visitation rates----
str(plant_pollinator_model_fit)
str(plant_with_pollinators)
treat <-levels(as.factor(plant_pollinator_model_fit$treatment))
treat <- treat[!treat %in% treatment.no.link]
treat <- treat[!treat %in% treatment.with.link]
# pol_visit_rate is create in HOI_1_Fusion_Dataframe
pol_visit_rate <- pol_visit_rate_0  %>%
    mutate(network=case_when(treatment  %in% c(treatment.no.link,"F","F_low","FF")  ~  "Link physically prevented",
                             treatment %in% c(treatment.with.link,"F","F_low","FF","NO_P") ~  "Nested structure" )) %>%
    mutate(focal=case_when(focal == "Raphanus" ~  "Radish",
                           focal == "Vicia"  ~  "Field bean" ,
                           focal== "Tomato" ~  "Tomato")) %>%
    mutate(species =case_when(species== "Bombus_terrestris" ~  "Bumblebee",
                              species == "Osmia_bicornis"  ~  "Bee",
                              species == "Lucilia_sericata"  ~  "Fly"
    ))

pol_visit_rate <-  arrange(transform(pol_visit_rate ,
                                     species=factor(species,level=c("Bumblebee","Bee" ,"Fly"))),
                           species)
pol_visit_rate <-  arrange(transform(pol_visit_rate ,focal=factor(focal,level=c("Radish","Field bean","Tomato"))),
                           focal)


write.table( pol_visit_rate,
           file.path("HOIs_Lynxc/results",
                     "pol_visit_rate.csv"), 
           na = "NA", append = F,
           col.names =TRUE)

Visitation_rate <- ggplot(pol_visit_rate , aes(y=value, x=focal, fill=as.factor(species))) + 
    #geom_bar(stat="identity",alpha= 0.8,position =position_dodge(),width=0.8) +  
    geom_boxplot() +
    facet_wrap(network~.,ncol=1,nrow=2) +   
    theme_bw() +  ylab("Number of visits") +  xlab("Focal") + 
    scale_fill_d3(name="",
                  labels=c("Bumblebee", "Mason bee", "Green bottle fly"))+
    theme(axis.text.x = element_text(face = "italic"),
          legend.text = element_text(face = "italic"),
          strip.background=element_rect(fill="white", color="white"),
          strip.text = element_text(size=10))

ggsave("HOIs_Lynxc/results/Visitation_rate.pdf",
       Visitation_rate, height = 5.27, width = 8.69, units = "in")

################################################################################
# F. Procrustes analysis
################################################################################
source('HOIs_Lynxc/HOI_6_Procrustes.R')
