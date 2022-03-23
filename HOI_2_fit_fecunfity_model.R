
#### Fit_fecundity_model, modified from Mayfield and Stouffer, 2018
# the model 'negbin-mm' used the factor year as random effect 
# gamma
#pollinator
################################################################################################
# Description of the function
##############################################################################################

fit.fecundity.model <- function(
    data,
    type=c("negbin", "negbin-mm", "poisson", "linear","inverse"),
    ...
){
    #  ---- Preliminary operations ----  
    # make sure the user has provided an eligible model
    match.arg(type)
    

    # drop extraneous levels that could complicate the model-fitting code
    data$year <- as.factor(data$year)
    data$year <- droplevels(data$year)
    not.numerical.col <- c("focal","seeds","treatment", "year")
    
    # lets figure out who the potential competitors are
    competitors <- colnames(data)[!colnames(data) %in% not.numerical.col]
    
    # remove species that are never observed to co-occur
    for(sp in names(which(colSums(data[,competitors, drop=FALSE],na.rm = T)==0))){
        data[,sp] <- NULL
    }
    
    # lets figure out the observed competitors 
    competitors <- colnames(data)[!colnames(data) %in% not.numerical.col]
    
    ##########################################################################################################
    # Description of the different coefficient based on the ecological network
    ##########################################################################################################
    # A ---- For pairwise competition----
    # all combinaison of plant-plant is possible 
    all.alphas <- competitors[! competitors %in% pollinator.id]
    
    # B ---- For mutualistic effect----
    # Combinaison of all plant-pollinator possible 
    all.gama <- competitors[!competitors %in% plant.id]
    all.gama.typeII <- unlist(lapply(all.gama, function(x){paste0("I(1/(1 + ",x,"))")}))
    
    # C ---- HOIs plant-pollinator-pollinator ----
    if(length(all.gama)>1){
        interbetas.polinator.on.mutu <- combn(all.gama,2)
        interbetas.polinator.on.mutu.2 <- apply(interbetas.polinator.on.mutu,2,paste,collapse=":")
        
        # combine all betas together into a single variable
        HOI_poll_poll <- c(interbetas.polinator.on.mutu.2)
    }else{HOI_poll_poll <-  c() } # c() need to be replaced by intrabetas.polinator.on.mutu if soft HOI included
    
    # D ---- HOIs plant-pollinator-plant ----
    ### fit.betas.polinator.on.comp ###
    HOI_poll_plants <- as.vector(outer(all.gama,all.alpha, paste, sep=":"))
    imp.interaction <- c("Lucilia_sericata:Tomato","Osmia_bicornis:Tomato","Lucilia_sericata:Vicia")
    imp.interaction.no.link <- c("Bombus_terrestris:Raphanus")
    HOI_poll_plants <- HOI_poll_plants[!HOI_poll_plants %in% imp.interaction]
    if (network == "no_link"){
        HOI_poll_plants<- HOI_poll_plants[!HOI_poll_plants%in% imp.interaction.no.link ]
    }
    # D ---- HOIs plant-plant-plant ----
    
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
    # Model with only the direct interactions
    ##########################################################################################################
    
    # statistically select for a linear or type II response from pollinator 
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
    
    result.dredge.pol <- as.data.frame(coef(Models.pol)) ##Final model selection table
    
    for (n.poll  in 1:length(all.gama)){
        result.dredge.pol <- result.dredge.pol[which((is.na(result.dredge.pol[,all.gama[n.poll]])  &
                                                         !is.na(result.dredge.pol[,all.gama.typeII[n.poll]])) 
                                                    | (is.na(result.dredge.pol[,all.gama.typeII[n.poll]]) &
                                                           !is.na(result.dredge.pol[,all.gama[n.poll]]))),]
    }
    
    result.dredge.pol <- result.dredge.pol[1,]
    result.dredge.pol <- result.dredge.pol[,!is.na(result.dredge.pol)]
    
    result.dredge.pol.names <- names(result.dredge.pol)[which(names(result.dredge.pol) %in% c(all.gama, all.gama.typeII))]
    
    ##########################################################################################################
    # fit the model
    # note that the call depends on the model form hence the many ifs (and they are not nested for clarity and because they don't need to be)
    ########################################################################################################## ##################
    model.formula <- paste(base.model.formula,
                           paste0(all.alpha, collapse=" + "),
                           paste0(result.dredge.pol.names, collapse=" + "),
                           #paste0(HOI_poll_poll, collapse=" + "),
                           paste0(HOI_plants_plants, collapse=" + "),
                           paste0(HOI_poll_plants, collapse=" + "),
                           sep=" + ")
    
    # fit the negative binomial model as described in the main text
    if(type=='negbin'){
        m <- glm.nb(
            as.formula(model.formula),
            data=data,
            na.action=na.fail,
            control=glm.control(maxit=1000),
            ...
        )
    }
    
    # fit a poisson form which is equivalent to the negative binomial without the dispersion parameter
    if(type=='poisson'){
        m <- glm(
            as.formula(model.formula),
            data=data, 
            family='poisson',
            control=glm.control(maxit=1000),
            ...
        )
    }


    
    
    # fit the linear form
    if(type=="linear"){
        m <- glm(
            as.formula(model.formula),
            data=data,
            control=glm.control(maxit=1000),
            ...
        )
    }
    
    # note that the mixed model and inverse forms are more finicky and harder to get to converge
    # the methods used to fit them also require a non-singular model matrix which is an issue for the beta model (that has many coefficients that cannot be inferred)
    
    # if(type=="negbin-mm" || type=="inverse"){
    # identify linearly independent predictors in this model
    # mm <- model.matrix(as.formula(model.formula), data)
    
    # separate out the terms that don't have to do with per capita effects
    # mm.base <- colnames(mm)[which(!colnames(mm) %in% all.alphas & !colnames(mm) %in% all.betas)]
    
    # figure out which columns (i.e. coefficients) to remove to eliminate rank deficiencies using qr decomposition
    # qrmm <- qr(mm)
    
    # separate out the names of these predictors
    #reduced.predictors <- colnames(mm)[qrmm$pivot[seq(qrmm$rank)]]
    
    # remove anything that is included in the base model so that the model formula command below makes sense
    # reduced.predictors <- reduced.predictors[which(!reduced.predictors %in% mm.base)]
    
    # separate these into alphas and betas
    #reduced.alphas <- reduced.predictors[which(!reduced.predictors %in% all.betas)]
    #reduced.betas <- reduced.predictors[which(reduced.predictors %in% all.betas)]
    
    # put them back together in a sorted order (so that coefficients in the the model summaries resemble what R would normally produce)
    # reduced.predictors <- c(sort(reduced.alphas), sort(reduced.betas))
    
    # reconstruct a non-rank-deficient model formula
    #if(length(reduced.predictors) > 0){
    #    model.formula <- paste0(base.model.formula," + ", paste0(reduced.predictors, collapse=" + "))
    #}else{
    #    model.formula <- base.model.formula
    #}
    #}
    # fit the negative binomial mixed model (with year level random effect)
    if(type=='negbin-mm'){
        m <- glmmadmb(
            as.formula(model.formula),
            random = as.formula(~ (1/"random.variable")),
            data=data,
            family='nbinom',
            ...
        )
    }
    
    # fit the inverse form
    if(type=='inverse'){
        m<-glm2(
            as.formula(model.formula),
            data=data,
            family=gaussian(link="inverse"),
            control=glm.control(maxit=1000),
            ...
        )
    }
    ##########################################################################################################
    # Model with only one type of response for pollinator 
    ##########################################################################################################
    model.formula.type.I <- paste(base.model.formula,
                           paste0(all.alpha, collapse=" + "),
                           paste0(all.gama, collapse=" + "),
                           #paste0(HOI_poll_poll, collapse=" + "),
                           paste0(HOI_plants_plants, collapse=" + "),
                           paste0(HOI_poll_plants, collapse=" + "),
                           sep=" + ")
  
        
     model.type.I <-   glm.nb(
        as.formula(model.formula.type.I),
        data=data,
        na.action=na.fail,
        control=glm.control(maxit=1000),
        ...
    )
     model.formula.type.II <- paste(base.model.formula,
                                   paste0(all.alpha, collapse=" + "),
                                   paste0(all.gama.typeII, collapse=" + "),
                                   #paste0(HOI_poll_poll, collapse=" + "),
                                   paste0(HOI_plants_plants, collapse=" + "),
                                   paste0(HOI_poll_plants, collapse=" + "),
                                   sep=" + ")
    model.type.II <- glm.nb(
        as.formula(model.formula.type.II),
        data=data,
        na.action=na.fail,
        control=glm.control(maxit=1000),
        ...
    )
    
    ##########################################################################################################
    # return data.frame and model
    ########################################################################################################## ##################
    
    order_terms <- c("Intercept","year2017",all.alphas,result.dredge.pol.names,
                     "Raphanus:Vicia", "Raphanus:Tomato","Tomato:Vicia",
                     HOI_poll_plants)
                   
    coef <- setNames(coefficients(m), order_terms)
    
return(list(model=m, coeff= coef,
            model.type.I =  model.type.I ,
            model.type.II = model.type.II ))

}


