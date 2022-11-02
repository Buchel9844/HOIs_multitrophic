# Code of Saavedra 2020, 
# Computes the fraction of LV simulations in which
# each individual species and each combination of species 
# survives and plots the results for a theoretical matrix
if(!require(deSolve)) { install.packages("deSolve"); library(deSolve) }
if(!require(diffeqr)) { install.packages("diffeqr"); library(diffeqr) }
if(!require(JuliaCall)) { install.packages("JuliaCall"); library(JuliaCall) }
library(MASS)
# 1°) Run the 3 function
# 2°) Simulation 
# 3°) plot
#############################
# functions
#############################
# Solves the system of ordinary differential equations
# given by the generalized Lotka-Volterra dynamics and
# returns the state variables over time


lotka_volterra <- function(N,r, K, Imat, times, formalism) {
    if (formalism == "r") {
        # list of parameters
        pars <- list(r=r,Imat=Imat)
        # function that returns rate of change
        model <- function(t, N, pars) {
            dN_dt <- N * (pars$r + c(pars$Imat %*% N))
            return(list(dN_dt))
        }
        # numerical integration
        out <- ode(y = N, times = times, func = model, parms = pars)
    }
    if (formalism == "K") {
        pars <- list(r = r, K = K,Imat = Imat)
        model <- function(t, N, pars) {
            dN_dt <- N * (pars$r / pars$K) * (pars$K + c(pars$Imat %*% N))
            return(list(dN_dt))
        }
        out <- ode(y = N, times = times, func = model, parms = pars)
    }
    if (formalism == "r_typeII") {
        pars <- list(r = r, Imat = Imat)
        model <- function(t, N, pars) {
            dN_dt <- N * (pars$r + c(pars$Imat %*% diag(1 / (1 + N)) %*% N))
            return(list(dN_dt))
        }
        out <- ode(y = N, times = times, func = model, parms = pars)
    }
    if (formalism == "r_stochastic") {
        pars <- list(r = r, Imat = Imat)
        # defining deterministic part
        f <- function(u, p, t) {
            deterministic <- u * (p$r + c(p$Imat %*% u))
            return(deterministic)
        }
        # defining stochastic part
        g <- function(u, p, t) {
            s <- rep(1 / sqrt(length(p$r)), length(p$r))
            stochastic <- s * u * (p$r - c(p$Imat %*% u))
            return(stochastic)
        }
        # integration time steps
        time_step <- times[2] - times[1]
        # numerical integration
        sol <- sde.solve(f = f, g = g, u0 = N, tspan = range(times), p = pars, saveat = time_step)
        out <- as.data.frame(cbind(sol$t, sol$u))
        names(out) <- paste("time", 1:length(pars$r))
    }
    return(out)
}


# Samples m vectors randomly on a n-dimensional simplex
simplex_sampling <- function(m, n) {
    r <- list()
    for (j in 1:m) {
        dist <- c(sort(runif(n-1, 0, 1)), 1)
        r[[j]] <- c(dist[1], diff(dist))
    }
    return(r)
}

# Samples m matrix of nxn randomly 
matrix_sampling <- function(A,B,n) {
    #Ilist <- list()
    #for (sp in 1:m) {
        Imat <- matrix(ncol=n,nrow=n)
        for (spi in 1:n) {
            for (spj in 1:n) {
                #extrem <- sort(c(A[spi, spj] , A[spi,spj] + B[spi,spj]))
  sample.mat <- rnorm(1e6,A[spi,spj] , abs(B[spi,spj]/3))
  p1 <- sample.mat[sample.mat < 5 & sample.mat > - 5]
  Imat[spi,spj] <- sample(p1, 1, replace = T, prob = pnorm(p1))
  
          if(spi==spj &  Imat[spi,spj] >0){
          p1 <- sample.mat[sample.mat<=0 & sample.mat > - 5]
          if(length(p1)==0){print("intra positive, redo sampling")   
            }
          Imat[spi,spj] <- sample(p1, 1, replace = T, prob = pnorm(p1))
          }
            }
        }
    return(Imat)
}

# Solves the system of ordinary differential equations
# given by the generalized Lotka-Volterra dynamics and
# returns the surviving species at equilibrium
lv_pruning <- function(N,r,Imat, K, times, formalism, extinct_tol = 0.00000001) {
    out <- lotka_volterra(N,r, K, Imat, times, formalism)
    eq <- out[nrow(out), -1]
    which_surv_sp <- as.numeric(which(eq > extinct_tol))
    return(paste(which_surv_sp, collapse = "_"))
}

#############################
# Simulation 
#############################
# setting simulation parameters ------------------------------
# number of initial conditions
init_cond <- 100
# sampling initial conditions randomly - Initial growth rate
n_sp <- 3
set.seed(1)
N <- replicate(init_cond, runif(n_sp, 0, 1), simplify = FALSE)
# number of simulations
n_sim_r <- 100
n_sim_alpha <- 10
# time steps
times <- seq(0, 200, 0.01)
# choose formalism
formalism <- "r" # other options are K (for K-formalism), "r_typeII, r_stochastic 

# loop to run the simulation for the different model and structure ------------------------------
# takes time, even for init_cond <- 10
noStd = T

for(HOIs in  c("with_HOIs","no_HOIs","full_model")){
    for(network in  c("with_link", "no_link")){

        print(paste0(HOIs,network))
        
        A <- as.matrix(Coefficients_alpha_r[[paste(network,HOIs,sep='_')]]$direct_int)
        if(HOIs == "no_HOIs"){
            B <- as.matrix(Coefficients_alpha_r[[paste(network,HOIs,sep='_')]]$direct_int_sd)
        }else{
        B <- as.matrix(Coefficients_alpha_r[[paste(network,HOIs,sep='_')]]$matrix_HOIs)
        }

        # maximum values : A+B
        
        
        #A <- as.matrix(Coefficients_alpha_r$Model.HOI.plant.on.plant.and.poll.on.mutualism_with_link$alpha)
        
        
        #A <- jitter(A,factor=5, amount=NULL)
        # add jitter to the smallest value 
        #A[which(A>0)] <- 0 
        # extracting number of species
        n_sp <- nrow(A)
        rownames(A) <- paste("sp", 1:n_sp, sep = "")
        colnames(A) <- paste("sp", 1:n_sp, sep = "")
        
        # performing simulations ------------------------------
        df <- data.frame()
        for (i in 1:init_cond) {
            print(i)
            # sample r vectors on the simplex
            r <- simplex_sampling(n_sim_r, n_sp)
            # sample m matrix of nxn randomly with a minimum of 2 inter-specific negative interactions
            Ilist <- list()
            counter <-1
            for(m in 1:100){
              Imat <- matrix_sampling(A,B, n_sp)
              if(sum(as.vector(sign(  Imat))) <= -1){
                Ilist[[counter]] <- Imat 
                counter <- counter + 1
              }
              if(counter > 20) break
            }
            #Imat <- matrix_sampling(A,B,50, n_sp)
            #Imat <- sample(Imat[sapply(Imat, function(x) all( sum(as.vector(sign(x))) <= -1, na.rm = TRUE))],11)
            # run simulations and extract surviving species
          surv_sp <- list()
          counter <- 1
            for (n in 1:20){
                # to remove symmetric matrix that have a_ij and a_ji positive, otherwise species i and j 
                if((length(Ilist[[n]][which(Ilist[[n]] > 0)]) >= 2 & isSymmetric(sign(Ilist[[n]])) )) next
              #print(n)
              
              tryCatch(surv_sp[[counter]] <- mapply(lv_pruning, r, MoreArgs = list( Imat=as.matrix( Ilist[[n]]),
                                                               N = N[[i]], times = times,
                                                               formalism = formalism,
                                                               extinct_tol = 0.0001),
                               SIMPLIFY = FALSE),
                        error = function(e) { skip_to_next <<- TRUE})
              
              counter <- counter + 1
              
              if(counter > n_sim_alpha)break
              
              # the 11th to 20th matrix are the back up matrixes in case the randomly choosen matrix leads 
              # to an mathematically unsolvabled ODE
            }
            
            # build data frame with results
            number_sim_real <- length(unlist(surv_sp))/n_sim_alpha
            r_df <- data.frame(matrix(rep(unlist(r),times=number_sim_real), nrow = length(r)*number_sim_real, byrow = TRUE))
            names(r_df) <- paste("r", 1:nrow(A), sep = "")
            N_df <- as.data.frame(matrix(rep(N[[i]], each = nrow(r_df)), ncol = nrow(A)))
            names(N_df) <- colnames(A)
            curr_df <- cbind(r_df, N_df)
            curr_df$init_cond <- i
            curr_df$surv_sp <- unlist(surv_sp)
            curr_df$perc <- sum(unlist(surv_sp) > 0)/ length(unlist(surv_sp))
            nrow(curr_df)
            curr_df <- curr_df[which(curr_df$surv_sp != ""),]
            # merge with larger data frame
            df <- rbind(df, curr_df)
        }
        
        
        # save results ------------------------------
        write.table(df, paste0("HOIs_Lynxc/results/predicting_surv_sp_",paste(network,HOIs,sep='_'), 
                             ".csv", sep = ""), 
                  row.names = FALSE)
    
    }
}

for(network in  c("with_link", "no_link")){
  HOIs= "no_HOIs"
  print(paste0(HOIs,network))
  
  A <- as.matrix(Coefficients_alpha_r[[paste(network,HOIs,sep='_')]]$direct_int)

  # maximum values : A+B
  
  
  #A <- as.matrix(Coefficients_alpha_r$Model.HOI.plant.on.plant.and.poll.on.mutualism_with_link$alpha)
  
  
  #A <- jitter(A,factor=5, amount=NULL)
  # add jitter to the smallest value 
  #A[which(A>0)] <- 0 
  # extracting number of species
  n_sp <- nrow(A)
  rownames(A) <- paste("sp", 1:n_sp, sep = "")
  colnames(A) <- paste("sp", 1:n_sp, sep = "")
  
  diag(A)[which( diag(A) >0)] <- -(abs(jitter(0)))
  # performing simulations ------------------------------
  df <- data.frame()
  for (i in 1:init_cond) {
    print(i)
    # sample r vectors on the simplex
    r <- simplex_sampling(n_sim_r, n_sp)
    
     # run simulations and extract surviving species
    surv_sp <- list()
      
      surv_sp <- mapply(lv_pruning, r, MoreArgs = list( Imat=A, N = N[[i]], times = times,
                                                                            formalism = formalism,
                                                                            extinct_tol = 0.0001),
                                            SIMPLIFY = FALSE)
 
    
    # build data frame with results
    number_sim_real <- length(unlist(surv_sp))/n_sim_alpha
    r_df <- data.frame(matrix(rep(unlist(r),times=number_sim_real), nrow = length(r)*number_sim_real, byrow = TRUE))
    names(r_df) <- paste("r", 1:nrow(A), sep = "")
    N_df <- as.data.frame(matrix(rep(N[[i]], each = nrow(r_df)), ncol = nrow(A)))
    names(N_df) <- colnames(A)
    curr_df <- cbind(r_df, N_df)
    curr_df$init_cond <- i
    curr_df$surv_sp <- unlist(surv_sp)
    curr_df$perc <- sum(unlist(surv_sp) > 0)/ length(unlist(surv_sp))
    nrow(curr_df)
    curr_df <- curr_df[which(curr_df$surv_sp != ""),]
    # merge with larger data frame
    df <- rbind(df, curr_df)
  }
  
  
  # save results ------------------------------
  
  str(df)
    write.table(df, 
                 file.path(paste0("HOIs_Lynxc/results/predicting_surv_sp_noStd_",network,HOIs, 
                           ".csv", sep = "")),
                 na = "NA", append = F,
                 col.names =TRUE)
}
view(df)
#############################
# Plot
#############################

# ggtern not compatibe with latest version of ggplot
# correct version are :  4.0.2 ggplot2: 3.3.2 ggtern: 3.3.0 - ignore the warning...
# the use of ggtern can be tricky, if showing recurrent warning, try to relaunch Rstudio. 
library(ggplot2)
library(ggtern)
# you need to relaod the packages of tidyverse separetely to avoid updtaing ggplot
library(tidyr)
library(plyr)
level.surv_sp <- c("Radish" ,"Field  bean","Tomato","Radish , Field  bean , Tomato",
                   "Radish , Field  bean","Radish , Tomato",
                   "Field  bean , Tomato" )
#pal_set1 <- c("#E41A1C", "#377EB8", "#FFFF33", "#984EA3", "#FF7F00", "#4DAF4A", "#A65628")
#pal_set1 <- c("#88CCEE", "#332288", "#6699CC", "#117733", "#CC6677", "#AA4499", "#882255")

summary.prob.sp<- NULL
predicting_surv_sp <- list()
predicting_surv_sp_colorblind <- list()
colorblind <- F # to make colorblind version
noStd <- T # to make the no hois without standard deviation version
for(network in  c("with_link", "no_link")){
    for(HOIs in  c("with_HOIs","no_HOIs","full_model")){
     # pal_set1 <- c("#88CCEE", "#332288", "#6699CC", "#117733", "#CC6677", "#AA4499", "#882255")
      pal_set1 <- c( "#00BFFF","#6699CC", "#79CDCD", "#FFD700", "#E066FF", "#8B008B", "#FF3E96")
      
       if(colorblind ==T){
      pal_set1<- pal_colorblind <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#000000")
      }
        print(paste0(HOIs,network))
        save_plots <- TRUE
        
        # reading results file ------------------------------
        head(   df_classic)
        df_classic <- read.csv(paste0("HOIs_Lynxc/results/predicting_surv_sp_", paste0(network,"_",HOIs), 
                                      ".csv"),sep=' ')
        if(noStd==T){
          if(!HOIs == "no_HOIs") next
          df_classic <- read.csv(paste0("HOIs_Lynxc/results/predicting_surv_sp_noStd_",network,HOIs, 
                                        ".csv"),sep=' ')
          network <- paste0("noStd_",network)
        }
        df_classic <-as.data.frame(df_classic)
        #head(df_classic)
      
        #sp_names <- paste("sp", 1:3, sep = "")
        sp_names <- c("Radish", "Field  bean", "Tomato")
        
        
        # classic (type I): building summary data frame for plotting ------------------------------
        # change combination names
        df_classic$surv_sp <- as.character(df_classic$surv_sp)
        df_classic$surv_sp[df_classic$surv_sp == "1"] <- sp_names[1]
        df_classic$surv_sp[df_classic$surv_sp == "2"] <- sp_names[2]
        df_classic$surv_sp[df_classic$surv_sp == "3"] <- sp_names[3]
        df_classic$surv_sp[df_classic$surv_sp == "1_2"] <- paste(sp_names[1],sp_names[2], sep=" , ")
        df_classic$surv_sp[df_classic$surv_sp == "1_3"] <- paste(sp_names[1],sp_names[3], sep=" , ")
        df_classic$surv_sp[df_classic$surv_sp == "2_3"] <- paste(sp_names[2],sp_names[3], sep=" , ")
        df_classic$surv_sp[df_classic$surv_sp == "1_2_3"] <- paste(sp_names[1],sp_names[2],sp_names[3], sep=" , ")
        df_classic$surv_sp <- factor(df_classic$surv_sp, levels = c(sp_names[1],sp_names[2],sp_names[3],
                                                                    paste(sp_names[1],sp_names[2],sp_names[3], sep=" , "),
                                                                    paste(sp_names[1],sp_names[2], sep=" , "),
                                                                    paste(sp_names[1],sp_names[3],sep=" , "),
                                                                    paste(sp_names[2],sp_names[3],sep=" , ")
                                                                   ))
        # build data frame with fraction of simulations with each combination
        init_cond <- unique(df_classic$init_cond)
        sp_counts <- tapply(df_classic$surv_sp, df_classic$init_cond, table)
        comb_df <- data.frame(matrix(unlist(sp_counts), nrow = length(sp_counts), byrow = TRUE))
        names(comb_df) <- names(sp_counts[[1]])
        comb_df <- comb_df / apply(comb_df, 1, sum)
        # build data frame for individual species
        sp1_df <- comb_df[ , grep(sp_names[1], names(comb_df))]
        sp1_sum_df <- data.frame(sp = sp_names[1], frac = apply(sp1_df, 1, sum), init_cond = init_cond)
        sp2_df <- comb_df[ , grep(sp_names[2], names(comb_df))]
        sp2_sum_df <- data.frame(sp = sp_names[2], frac = apply(sp2_df, 1, sum), init_cond = init_cond)
        sp3_df <- comb_df[ , grep(sp_names[3], names(comb_df))]
        sp3_sum_df <- data.frame(sp = sp_names[3], frac = apply(sp3_df, 1, sum), init_cond = init_cond)
        sp_sum_df_classic <- rbind(sp1_sum_df, sp2_sum_df, sp3_sum_df)
        # data frame with all species combinations
        comb_df <- gather(comb_df, "surviving_sp")
        comb_df$init_cond <- rep(init_cond, length(unique(comb_df$surviving_sp)))
        names(comb_df) <- c("surviving_sp", "frac_solutions", "init_cond")
        summ_df_classic <- ddply(comb_df, "surviving_sp", summarise,
                                 mean_frac = mean(frac_solutions),
                                 sd_frac = sd(frac_solutions))
        
        
        # build full data frame with all types of dynamics ------------------------------
        sp_sum_df <- sp_sum_df_classic
        
        write.csv(sp_sum_df, paste0("HOIs_Lynxc/results/probability_surv_sp_", paste(network,HOIs,sep='_'), 
                                    ".csv", sep = ""), 
                  row.names = FALSE)
        write.csv(df_classic , paste0("HOIs_Lynxc/results/probability_df_classic_", paste(network,HOIs,sep='_'), 
                                    ".csv", sep = ""), 
                  row.names = FALSE)
        summary.prob.sp.i <- data.frame(summary(df_classic$surv_sp)/nrow(df_classic))
        names(summary.prob.sp.i) <- paste(network,HOIs,sep=' ')
        summary.prob.sp <- bind_cols(summary.prob.sp,summary.prob.sp.i)
        # creating palette
        # using just one initial condition
        # single_init_cond <- 1
        # data frame for the single initial condition
        # df_classic_init_cond <- as.data.frame(subset(df_classic, init_cond == single_init_cond))
        if (sum(summ_df_classic$mean_frac ==0)){
          pal_set1  <- pal_set1[!level.surv_sp %in% 
                                  summ_df_classic$surviving_sp[which(summ_df_classic$mean_frac ==0)]]
        }
        fig_A <- ggtern(data = df_classic, 
                        mapping= aes(x=r2, 
                                     y=r3,
                                     z=r1,
                                     color = surv_sp)) +
            geom_point(size = 1,shape=21,alpha=0.8) + 
            scale_color_manual(values = pal_set1,
                               #limit=level.surv_sp, 
                               name = "Observed\nensemble") +
            guides(color = guide_legend(override.aes = list(stroke = 2,size=4,alpha=1),
                                        nrow=2,byrow=TRUE)) +
            labs(x = expression(paste("Field \nbean")), y = "Tomato", z = "Radish") +
            #theme_bw() +
            theme_noticks() +
            theme(text = element_text(size=12), 
                  panel.grid.major = element_line(color = "white"),  
                  panel.grid.minor = element_line(color = "white"), 
                  panel.border = element_rect(color = "white", fill = "white"),
                  #legend.position = "none",
                  axis.text=element_text(size=5,color='white'),
                  axis.text.y = element_text(size = 12),
                  plot.background = element_rect(color = "white", fill = "white"),
                  axis.title = element_text(size = 12),
                  axis.text.x = element_text(size = 12),
                  legend.key = element_blank(),
                  legend.position="bottom"
            ) 
        #annotate(plot(Interaction_network[[paste(model,network,sep="_")]]),
        #           xmin = 1, xmax=1, ymin=1, ymax=1, zmin=1, zmax=1)
        
        # transform the figure into a ggplotgrop in order to plot all the triangle together
        # note: here ggarange does not work bc we work with ggterm and not ggplot
        if(colorblind == T){
          predicting_surv_sp_colorblind[[paste(network,HOIs,sep="_")]] <- fig_A
          ggsave(paste0("HOIs_Lynxc/results/Colorblind_",network,HOIs,".png"),fig_A,
                 width = 5.8,height = 5.2, units = "in")
        }else{
        predicting_surv_sp[[paste(network,HOIs,sep="_")]] <- fig_A
        ggsave(paste0("HOIs_Lynxc/results/",network,HOIs,".png"),fig_A,
               width = 5.8,height = 5.2, units = "in")
        }
        
    
    }
}

write.csv(  summary.prob.sp , paste0("HOIs_Lynxc/results/summary.prob.sp", 
                              ".csv", sep = ""), 
          row.names = T)
plot(predicting_surv_sp[[paste(network,HOIs,sep="_")]])
View(df_classic)
# get the legend  to create a common legend 
legend.plot <- ggplot(data = df_classic, 
                      mapping= aes(x = r2, 
                                   y = r3,
                                   z = r1,
                                   color =surv_sp)) + geom_point() +    theme_bw()  +
    guides(color = guide_legend(override.aes = list(stroke = 2,shape=21, size=4,alpha=1),
                              nrow=2,byrow=TRUE)) +
    theme(legend.position = "bottom",
          plot.margin = margin(0, 0, 0, 0, "cm")) + 
    scale_color_manual(values = pal_set1, #limit=levels(df_classic$surv_sp),
                       name = "Observed\nensemble")

ggsave("HOIs_Lynxc/results/predicting_surv_sp_legend.png",
       plot=plot(get_legend(legend.plot)), 
       width = 7.53, height =3.68, units = "in")

legend.plot.colorblind <- ggplot(data = df_classic, 
                      mapping= aes(x = r2, 
                                   y = r3,
                                   z = r1,
                                   color =surv_sp)) + geom_point() +    theme_bw()  +
  guides(color = guide_legend(override.aes = list(stroke = 2,shape=21,size=4,alpha=1),
                              nrow=2,byrow=TRUE)) +
  theme(legend.position = "bottom",
        plot.margin = margin(0, 0, 0, 0, "cm")) + 
  scale_color_manual(values = pal_colorblind, #limit=levels(df_classic$surv_sp),
                     name = "Observed\nensemble")



ggsave("HOIs_Lynxc/results/predicting_surv_sp_legend_colorblind.png",
       plot=plot(get_legend(legend.plot.colorblind )), units = "in")


library(ggpubr)
ggsave("HOIs_Lynxc/results/predicting_surv_sp.png",
       plot=ggarrange(plotlist=predicting_surv_sp, ncol=3,row=3,
                      common.legend = T,legend="bottom"), units = "in")

