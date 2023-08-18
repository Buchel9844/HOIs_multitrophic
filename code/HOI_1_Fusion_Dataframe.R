
################################################################################
# Import data
#################################################################################
#----plant seed and neighboors abundance----
plant <- read.csv("data/used_data2017.csv")
plant <- subset(plant,
                  select=c("link" ,"treatment","pot","focal",
                           "background_T","background_R","background_H",
                           "fruit_number","mean_weigthxseed","mean_seedsxfruit"
                           ,"total_seeds","total_seed_weigth"))

treatment <-  levels(as.factor(plant$treatment))

plant$focal[plant$focal =="T" ] <- "Tomato"
plant$focal[plant$focal =="H" ] <- "Vicia"
plant$focal[plant$focal =="R" ] <- "Raphanus"
plant$year <- "2017"
pollinator.id <- c("Bombus_terrestris","Lucilia_sericata","Osmia_bicornis")
for ( poll in pollinator.id){
    plant[,poll] <- 0
}


# Check 1
# str(plant) --> 738 obs. of  16 variables

#---- pollinator abundance----
plant_with_pollinators <- read.csv("data/pollinators.csv")
plant_with_pollinators <- data.frame(plant_with_pollinators[102:nrow(plant_with_pollinators),])

# separate the information in the column "measure"
plant_with_pollinators <- separate(plant_with_pollinators,measure,sep="_",
                        into = c("Number_visits", "visits", "focal", "time")) 
# only keep the interesting columns
plant_with_pollinators <- subset(plant_with_pollinators, 
                      select=c("year","species","focal","treatment"
                               ,"value"))
# make the column uniforme (same levels)
plant_with_pollinators$focal[plant_with_pollinators$focal =="Tomato" | 
                                 plant_with_pollinators$focal =="tomato" ] <- "Tomato"
plant_with_pollinators$focal[plant_with_pollinators$focal =="vicia" ] <- "Vicia"
plant_with_pollinators$focal[plant_with_pollinators$focal =="raphanus" ] <- "Raphanus"
pol_visit_rate_0 <- plant_with_pollinators

# make one columns that enables the regroupement by treatment 
plant_with_pollinators  <- unite(plant_with_pollinators, species,treatment,
                        col = "poll_treatment",sep = "-")
plant_with_pollinators  <- unite(plant_with_pollinators,poll_treatment,focal,
                      col = "poll_treatment_focal",sep = "-") 

# mean of the replicate by treatment
plant_with_pollinators <- plant_with_pollinators %>% 
                            group_by(poll_treatment_focal)%>%
                            dplyr::summarize(Mean = mean(value, na.rm=TRUE)) 

plant_with_pollinators$year <- "2017"
names(plant_with_pollinators) <- c("poll_treatment_focal","Mean","year")
# separation of the informations
plant_with_pollinators <- separate(plant_with_pollinators, poll_treatment_focal, sep="-",
             into = c("pollinator", "treatment", "focal"))

# these are the vectors which will be used for species identity throught the entire scripts
# check is there are correct and correspond to the columns' names 
pollinator.id <- levels(as.factor(plant_with_pollinators$pollinator))
plant.id <- levels(as.factor(plant_with_pollinators$focal))
plant.id <- c("Raphanus","Vicia","Tomato")
# Check 2
# str(plant_with_pollinators) --> 45 obs. of  5 variables

#----plant seed and neighboors abundance wihtout pollinators----
plant_no_pollinator <- read.csv("data/2016_no_pollinators_for_2017.csv")

plant_no_pollinator <- subset(plant_no_pollinator,
                              select=c("link" ,"treatment","pot","focal",
                                       "background_T","background_R","background_H",
                                       "fruit_number","mean_weigthxseed","mean_seedsxfruit"
                                       ,"total_seeds","total_seed_weigth"))

for ( poll in pollinator.id){
    plant_no_pollinator[,poll] <- NA
}

plant_no_pollinator$focal[plant_no_pollinator$focal =="T" ] <- "Tomato"
plant_no_pollinator$focal[plant_no_pollinator$focal =="H" ] <- "Vicia"
plant_no_pollinator$focal[ plant_no_pollinator$focal =="R" ] <- "Raphanus"
plant_no_pollinator$year <- "2016"

# Check 3
# str(plant_no_pollinator) --> 168 obs. of  16 variables

################################################################################
# Combinaison of plant and pollinator
#################################################################################

for ( spcs in plant.id){
    for ( treat in treatment){
        for ( poll in pollinator.id){
            if (length(plant_with_pollinators$Mean[plant_with_pollinators$pollinator == poll &
                                                   plant_with_pollinators$focal ==  spcs &
                                                   plant_with_pollinators$treatment == treat])==0) next
            else{
            plant[plant$focal ==  spcs  &
                      plant$treatment == treat, poll] <- plant_with_pollinators$Mean[plant_with_pollinators$pollinator == poll &
                                                                                         plant_with_pollinators$focal ==  spcs &
                                                                                         plant_with_pollinators$treatment == treat]
            }
            
        }
    }
}

# Check 3
# str(plant) --> 738 obs. of  16 variables


################################################################################
# Combinaison of plant and plant without pollinators
#################################################################################

plant_pollinator <- bind_rows(plant, plant_no_pollinator)

# Check 4
# str(plant_pollinator) --> 906 obs. of  16 variables

# "total_seeds" is the fitness approximation
plant_pollinator_model_fit <- subset(plant_pollinator, select=c("year","treatment","total_seeds","focal",
                                                                "background_T","background_R","background_H",
                                                                "Bombus_terrestris","Lucilia_sericata","Osmia_bicornis"))

names(plant_pollinator_model_fit) <- c("year","treatment", "seeds", "focal", plant.id, pollinator.id)

#plant_pollinator_model_fit$seeds <- ceiling(plant_pollinator_model_fit$seeds)
plant_pollinator_model_fit <- replace_na(plant_pollinator_model_fit,
                                        replace = list(Bombus_terrestris = 0,
                                                       Lucilia_sericata = 0,
                                                       Osmia_bicornis = 0))

# Check 5
# plant_pollinator_model_fit[which(is.na(plant_pollinator_model_fit) == T),] --> A tibble: 0 x 10

# we need to separate the treatment with and without link (see main text)
plant_pollinator_model_fit$treatment <- as.factor(plant_pollinator_model_fit$treatment)
levels(as.factor(plant_pollinator_model_fit$treatment))
treatment.no.link <- c("B_x","BB_x","BF_x","OB_x")
treatment.with.link <- c("B","B_low","BB","BF","OB","OF","OO","O","O_low")
# the treatments ("F","F_low","FF","NO_P") should not involve different result between nested or half nester structure
plant_pollinator_model_fit$treatment <- as.character(plant_pollinator_model_fit$treatment)
# creation of a list of the data.frame for each focal species with the full("with_link") and partial network("no_link")
fecundity.data <- list()

#plant_pollinator_model_fit$seeds <- round(plant_pollinator_model_fit$seeds)
maxplant <- max(plant_pollinator_model_fit$Raphanus,
                plant_pollinator_model_fit$Vicia,
                plant_pollinator_model_fit$Tomato)
maxpol <-max(plant_pollinator_model_fit$Bombus_terrestris,
    plant_pollinator_model_fit$Lucilia_sericata,
    plant_pollinator_model_fit$Osmia_bicornis)

maxpol.log <- max(log(plant_pollinator_model_fit$Bombus_terrestris),
              log(plant_pollinator_model_fit$Lucilia_sericata),
              log(plant_pollinator_model_fit$Osmia_bicornis))

VR.plot <- list()
for(plant.spc in plant.id){
 df.of.focal <-  as.data.frame(subset(plant_pollinator_model_fit,
                         plant_pollinator_model_fit$focal == plant.spc))
 #df.of.focal$Raphanus <- df.of.focal$Raphanus/maxplant
 #df.of.focal$Vicia <- df.of.focal$Vicia/ maxplant
 #df.of.focal$Tomato <- df.of.focal$Tomato/ maxplant
 df.of.focal.test <-  df.of.focal 
 df.of.focal.test$Bombus_terrestris.log <- (log(df.of.focal$Bombus_terrestris)/ maxpol.log)*maxplant
 df.of.focal.test$Lucilia_sericata.log <- (log(df.of.focal$Lucilia_sericata)/ maxpol.log)*maxplant
 df.of.focal.test$Osmia_bicornis.log <- (log(df.of.focal$Osmia_bicornis)/maxpol.log)*maxplant
 
 df.of.focal.test$Bombus_terrestris.log[which( df.of.focal.test$Bombus_terrestris.log < -1)] <- -1
 df.of.focal.test$Bombus_terrestris.log <-  df.of.focal.test$Bombus_terrestris.log +1
 df.of.focal.test$Lucilia_sericata.log[which( df.of.focal.test$Lucilia_sericata.log < -1)] <- -1
 df.of.focal.test$Lucilia_sericata.log <-   df.of.focal.test$Lucilia_sericata.log + 1
 df.of.focal.test$Osmia_bicornis.log[which( df.of.focal.test$Osmia_bicornis.log < -1)] <- -1
 df.of.focal.test$Osmia_bicornis.log <-   df.of.focal.test$Osmia_bicornis.log +1
 
 
 df.of.focal.test$Bombus_terrestris.std <- (df.of.focal$Bombus_terrestris/ maxpol)*maxplant
 df.of.focal.test$Lucilia_sericata.std <- (df.of.focal$Lucilia_sericata/ maxpol)*maxplant
 df.of.focal.test$Osmia_bicornis.std <- (df.of.focal$Osmia_bicornis/maxpol)*maxplant
 
 
 VR.plot[[plant.spc]] <- ggplot(df.of.focal.test) + 
   geom_smooth(aes(y=seeds, x=(Bombus_terrestris+
                                 Lucilia_sericata+
                                 Osmia_bicornis),
                   colour= "1 - Original distribution V")) +
   
   geom_smooth(aes(y=seeds, x=(Bombus_terrestris.std+
                                Lucilia_sericata.std+
                                Osmia_bicornis.std),
                   colour="2- {V divided by max(V)}* max(plant abundance)")) + 
   #xlim(c(0,5)) + 
   #ylim(c(0,100)) + 
   geom_smooth(aes(y=seeds, x=(Bombus_terrestris.log+
                                 Lucilia_sericata.log+
                                 Osmia_bicornis.log),
                   colour="3- {log(V) divided by max(log(V))}* max(plant abundance)")) +
   scale_colour_manual(name="standardisation", values=c("blue", "orange","black")) +
   theme_bw() + ggtitle(plant.spc) +
   xlab("Visitations rate")
 
 
 df.of.focal$Bombus_terrestris <-  df.of.focal.test$Bombus_terrestris.log
 df.of.focal$Lucilia_sericata <-  df.of.focal.test$Lucilia_sericata.log
 df.of.focal$Osmia_bicornis <-  df.of.focal.test$Osmia_bicornis.log
 
 fecundity.data[[paste(plant.spc,"with_link", sep="_")]] <- df.of.focal[!df.of.focal$treatment %in% treatment.no.link, ]
 fecundity.data[[paste(plant.spc,"no_link", sep="_")]] <- df.of.focal[!df.of.focal$treatment %in% treatment.with.link, ]
 
}

save(fecundity.data,
     file = "results/fecundity_data.RData")


ggplot2::ggsave("results/VisitationRate_standardisation.pdf",
  ggarrange(plotlist = VR.plot, common.legend = T, legend="bottom"),
  width = 11.20,height = 5.56, units = "in")
# Check 6
# levels(as.factor(fecundity.data$Raphanus_no_link$treatment)) --> c("B_x","BB_x","BF_x","F","F_low","FF","NO_P","OB_x")
# OR length(levels(as.factor(fecundity.data$Raphanus_no_link$treatment))) -->  8
################################################################################
# Save new data frame
#################################################################################
View(plant_pollinator)
write.table(plant_pollinator,
          file.path("results/plant_pollinator_2016_2017.csv"), 
          na = "NA", append = F,
          col.names =TRUE)
fecundity.data$Raphanus_with_link


