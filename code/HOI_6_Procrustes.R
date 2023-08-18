# Simple Procruste Rotation between two sets of points
# DETAILS: 
# Here X and y are not independant, so it does not matter which matrix is X and which is Y. 
# Function protest is a permutational test of the significance of the procrustes
# result. It is based on the correlation from a symmetric Procrustes analysis. 
# This is why the order of X and Y makes no difference to the test of significance.
# Null hypothesis: 	The degree of concordance between two (or more) matrices
# is no greater than expected given random inter-matrix associations.
# exampl: if pvalue (significance here) = 0.08
# The p-value suggests that about 8 times in 100 you'd get the same test statistic (similar correlation) if you randomised one of the matrices.

# Procrustes statistic, m2 : The degree to which the superimposition was successful is determined
# smaller values of m12 indicate higher con- cordance between data sets 
# Significance of the m2 statistic is often tested by permutation (see Jackson, 1995), 
# whereby the row assignments in one matrix are randomly permuted a large number of times to create the null distribution
#  #the probability of rejection is assessed as (number of m12-rnd equal to or smaller than m12-obs+1)/(number of randomizations+1). The 1 in the numerator and the denominator represents the observed value for the statistic being evaluated, which is considered as a possible value of the randomized

# Test for diff betwenn no hoi and with hoi for both structure

############### Boot protest 

nsim <- 100
#set the number of simulation
    
boot.protest <- function(matrix.1,matrix.sd.1,
                         matrix.2,matrix.sd.2){
    sample.matrix.1 <- data.frame(Raphanus_Raphanus= 1:(nsim*10))
    sample.matrix.2 <- data.frame(Raphanus_Raphanus= 1:(nsim*10))
    for (i in plant.id){
        for (n in plant.id){
            sample.matrix.1[,paste(i,n,sep="_")] <- rnorm((nsim*10),mean=matrix.1[i,n],
                                                          sd=matrix.sd.1[i,n])
            sample.matrix.2[,paste(i,n,sep="_")] <- rnorm((nsim*10),mean=matrix.2[i,n],
                                                          sd=matrix.sd.2[i,n])
        }
    }
    scale <- c()
    for (row.i in sample(1:(nsim*10),nsim)){
        row.matrix.1 <- as.numeric(sample.matrix.1[row.i,])
        alpha.matrix.1 <- t(matrix(row.matrix.1, ncol=3, nrow=3))
        row.matrix.2 <- as.numeric(sample.matrix.2[row.i,])
        alpha.matrix.2 <- t(matrix(row.matrix.2, ncol=3, nrow=3))
            
        fit <- vegan::procrustes(alpha.matrix.1,alpha.matrix.2,
                                 scale = TRUE, symmetric = T, 
                                 permutation= how(nperm = 999))
        scale <- c(scale,fit$scale)
        }
    return(scale)
}


for(network in  c("with_link", "no_link")){
    for(HOIs0 in  c("no_HOIs", "with_HOIs","full_model")){
    for(HOIs in  c("no_HOIs", "with_HOIs","full_model")){
        model1 <- paste(network,HOIs0,sep="_")
        model2 <- paste(network,HOIs,sep="_")
        boot.df <- data.frame(network = rep(network,
                                              times=nsim),
                              model0= rep(HOIs0,
                                          times=nsim),
                              model1= rep(HOIs,
                                          times=nsim),
                              model.comb = rep(paste(HOIs0,HOIs,sep="_"),
                                            times=nsim))
        boot.df[,"Correlation"] <- boot.protest(Coefficients_alpha_r[[model1]]$alpha,
                                    Coefficients_alpha_r[[model1]]$sd,
                                    Coefficients_alpha_r[[model2]]$alpha,
                                    Coefficients_alpha_r[[model2]]$sd)
        if( network=="with_link" & HOIs0=="no_HOIs" & HOIs=="no_HOIs"){
            boot.model <- boot.df
        }else{  boot.model <- bind_rows(boot.model,boot.df)}
       
    }
    }
}

boot.model <- arrange(transform(boot.model,
                                 network =factor(network ,level=c("no_link","with_link"))),
                       network )

boot.model <- arrange(transform(boot.model,
                                 model0=factor(model0,level=c("no_HOIs","with_HOIs","full_model"))),
                       model0)
boot.model <- arrange(transform(boot.model,
                                model1=factor(model1,level=c("no_HOIs","with_HOIs","full_model"))),
                      model1)

boot.model <- boot.model %>%
    mutate(network =case_when(network== "no_link" ~  "Link physically prevented",
                          network == "with_link"  ~  "Nested structure"
)) %>%
    mutate(model0 =case_when(model0 == "no_HOIs" ~  "Pairwise interactions",
                            model0 == "with_HOIs"  ~  "Selected HOIs",
                            model0 == "full_model" ~  "All HOIs"
    ))%>%
    mutate(model1 =case_when(model1 == "no_HOIs" ~  "Pairwise interactions",
                             model1 == "with_HOIs"  ~  "Selected HOIs",
                             model1 == "full_model" ~  "All HOIs"
    ))

write.table( boot.model,
              "HOIs_Lynxc/results/boot.model.csv", 
              na = "NA", append = F)
    

boot.model <- arrange(transform(boot.model,
                                network =factor(network ,level=c("Link physically prevented",
                                                                 "Nested structure"))),
                      network )

boot.model <- arrange(transform(boot.model,
                                model0=factor(model0,level=c("Pairwise interactions",
                                                             "Selected HOIs",
                                                             "All HOIs"))),
                      model0)
boot.model <- arrange(transform(boot.model,
                                model1=factor(model1,level=c("Pairwise interactions",
                                                             "Selected HOIs",
                                                             "All HOIs"))),
                      model1)

density.procrustes.plot <- ggplot(boot.model, aes(x=Correlation)) +
    geom_density(alpha = 0.4, aes(color = model1), kernel="gaussian") +
    theme_bw()  + scale_color_npg(name= "") +
    theme(legend.position = "bottom") +
    facet_wrap( network ~ model0,
                labeller = labeller(.multi_line = T)) + 
    theme(strip.background=element_rect(fill="white", color="white"),
          strip.text = element_text(size=10))+
    xlab("Correlation between two interaction matrices")




ggsave("HOIs_Lynxc/results/density.procrustes.plot.pdf",
           plot=density.procrustes.plot,
           height = 8.27, width = 11.69, units = "in")

