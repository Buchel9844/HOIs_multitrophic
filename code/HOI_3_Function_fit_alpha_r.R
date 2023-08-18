#R-code of "Will a small randomly-assembled community be feasible and stable?" by:
#Chuliang Song and Serguei Saavedra
#published in: Ecology


require(mvtnorm)

Omega <- function(alpha){
    n <- nrow(alpha)
    if (abs(det(alpha)) > 1e-6){
        Sigma <-solve(t(alpha) %*% alpha)
        d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
        out <- log10(d[1]) + n * log10(2)
    } else {out <- NA}
    return(out) 
}

######################################

norm_L2 <- function(alpha){
    D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
    out <- alpha %*% D
    return(out)
}

######################################

r_centroid <- function(alpha){
    r_c <- rowSums(alpha)
    r_c / sqrt((sum(r_c^2)))
    r_c <- t(t(r_c))
    return(r_c)
}

######################################

theta <- function(v1,v2){
    out <- acos(sum(v1*v2)/(sqrt(sum(v1^2))*sqrt(sum(v2^2))))*180/pi
    return(out)
}

######################################
#test if a system (alpha and r) is feasible (output 1 = feasible, 0 = not feasible)
test_feasibility <- function(alpha,r){
    out <- prod(solve(alpha, r) > 0)
    return(out)
}

######################################
#test which pairs in a system (alpha and r) are feasible (output 1 = feasible, 0 = not feasible)
test_feasibility_pairs <- function(alpha,r){
    n <- length(r)
    c <- combn(n,2)
    nc <- dim(c)[2]
    f <- rep(NA,nc)
    for (i in 1:nc){
        f[i] <- prod(solve(alpha[c[,i],c[,i]],r[c[,i]])>0)
    }
    out <- list(pairs = c, feasibility = f)
    return(out)
}

######################################
#compute the feasiblity domain, the feasibility domain of all pairs, and their overlap (Nrand = number of randomization)
compute_overlap <- function(alpha,Nrand){
    
    n <- dim(alpha)[1]
    
    counter_f <- 0
    counter_overlap <- 0
    counter_all <- 0
    
    for (i in 1:Nrand){
        
        r_rand <- abs(rnorm(n))  
        r_rand <- r_rand/sqrt(sum(r_rand^2))
        
        f1 <- test_feasibility(alpha,r_rand)  
        f2 <- test_feasibility_pairs(alpha,r_rand)$feasibility  
        
        counter_f <- counter_f + f1
        counter_all <- counter_all + prod(f2)
        counter_overlap <- counter_overlap + f1*prod(f2)
        
    }
    
    Omega <- counter_f/Nrand
    Omega_all <- counter_all/Nrand
    overlap <- counter_overlap/Nrand
    Differential <- Omega_all - Omega
    
    out <- list(Omega = Omega, Omega_all = Omega_all, overlap = overlap,  Differential =  Differential)
    return(out)
    
}
