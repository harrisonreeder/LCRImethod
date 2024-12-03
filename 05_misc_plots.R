library(nleqslv)
palette <- RColorBrewer::brewer.pal(n=7,name = "Set2")


####SIMPLIFIED GRAPHIC####
cairo_pdf(filename = "pi_beta0_settings_simplified.pdf",
          width = 3.5,height=3.5)
palette <- RColorBrewer::brewer.pal(n=4,name = "Set2")
plot(x=0,y=0,pch=NA,type="p",xlim=c(0,0.5),ylim=c(0,1),axes=FALSE,
     xlab="",ylab="")
axis(side=1,at=c(0,0.25,0.333,0.5),las=2,cex.axis=0.8)
axis(side=2,at=c(0,0.25,0.5,0.75,1),las=1,cex.axis=0.8)
title(#main="main title", sub="sub-title",
  xlab=expression(PASC~Rate~among~Infected~(pi)),
  ylab=expression(Symptom~Prevalence~among~non-PASC~(beta[0])),cex.lab=0.7)
x <- seq(0,0.5,0.001)
y <- 1- x/(1-x)
lines(x,y)
polygon(x=c(x,0,0),y=c(y,0,1),col="grey70",#density=20,angle=90,
        border = NA)
x <- c(0,1/3);y <- c(1/2,1/2)
lines(x=x,y=y,lwd=1.5,col="black",lty=3)
x <- c(1/3,1/3);y <- c(0,1/2)
lines(x=x,y=y,lwd=1.5,col="black",lty=3)
lines(x=c(0,0),y=c(0,1),lwd=0.75,col="black",lty=1)
lines(x=c(0.5,0.5),y=c(0,1),lwd=0.75,col="black",lty=1)
lines(x=c(0,0.5),y=c(1,1),lwd=0.75,col="black",lty=1)
lines(x=c(0,0.5),y=c(0,0),lwd=0.75,col="black",lty=1)
legend(x="topright",legend=c(expression(paste(1/beta[0]) < phi),
                             expression(paste(1/beta[0]) >= phi)),
       fill=c("grey70","white"),cex=0.7)
dev.off()





####addition of confounder ----

#function that takes in a variety of inputs as described in the paper,
#and (numerically) computes the corresponding odds ratio, with and without balancing weights

OR_func <- function(pr_A, #alpha
                    pr_Y_A1, #pi
                    pr_X_Y0, #beta0
                    rr_X_Y, #beta1
                    pr_Z, rr_A_Z, rr_X_Z){
  # browser()
  #derived quantities
  # #https://www.wolframalpha.com/input?i=solve+a+%3D+b+*+c+*+d+%2B+%281-b%29+*+c+for+c
  pr_A_Z0 <- pr_A / (pr_Z * (rr_A_Z-1) + 1)
  pr_A_Z1 <- pr_A_Z0 * rr_A_Z
  pr_Z_A0 <- (1-pr_A_Z1) * pr_Z / (1-pr_A)
  pr_Z_A1 <- pr_A_Z1 * pr_Z / pr_A
  stopifnot(0<pr_A_Z0 & pr_A_Z0<1)
  stopifnot(0<pr_A_Z1 & pr_A_Z1<1)
  stopifnot(pr_A == pr_Z * pr_A_Z1 + (1-pr_Z) * pr_A_Z0)
  stopifnot(pr_Z == pr_A * pr_Z_A1 + (1-pr_A) * pr_Z_A0)
  
  #derived quantities
  pr_Y_A0 <- pr_Y_A0Z1 <- pr_Y_A0Z0 <- 0
  pr_Y_A1Z1 <- pr_Y_A1Z0 <- pr_Y_A1
  pr_Y <- pr_A * pr_Y_A1 + (1-pr_A) * pr_Y_A0
  pr_Y_Z0 <- pr_A_Z0 * pr_Y_A1Z0 + (1-pr_A_Z0) * pr_Y_A0Z0
  pr_Y_Z1 <- pr_A_Z1 * pr_Y_A1Z1 + (1-pr_A_Z1) * pr_Y_A0Z1
  pr_Z1Y0 <- pr_Z * (1-pr_Y_Z1)
  pr_Z1Y1 <- pr_Z * pr_Y_Z1
  pr_Z_Y0 <- pr_Z1Y0 / (1-pr_Y)
  pr_Z_Y1 <- pr_Z1Y1 / pr_Y
  rr_Y_Z <- pr_Y_Z1 / pr_Y_Z0
  
  #numerically find the conditional risk ratios between X and each of Y and Z, 
  #conditional on the other
  score_func <- function(x){
    pr_X_Y0Z0 <- x[1]
    rr_X_Y_Z <- x[2]
    rr_X_Z_Y <- x[3]
    
    pr_X_Y1 <- pr_X_Y0 * rr_X_Y
    pr_X_Y1Z0 <- pr_X_Y0Z0 * rr_X_Y_Z
    pr_X_Y0Z1 <- pr_X_Y0Z0 * rr_X_Z_Y
    pr_X_Y1Z1 <- pr_X_Y0Z0 * rr_X_Z_Y * rr_X_Y_Z
    pr_X_Z0 <- (1 - pr_Y_Z0) * pr_X_Y0Z0 + pr_Y_Z0 * pr_X_Y1Z0
    pr_X_Z1 <- (1 - pr_Y_Z1) * pr_X_Y0Z1 + pr_Y_Z1 * pr_X_Y1Z1
    c(
      rr_X_Z - pr_X_Z1/pr_X_Z0,
      pr_X_Y0 - (pr_Z_Y0 * pr_X_Y0Z1 + (1-pr_Z_Y0) * pr_X_Y0Z0),
      pr_X_Y1 - (pr_Z_Y1 * pr_X_Y1Z1 + (1-pr_Z_Y1) * pr_X_Y1Z0)
    )
  }
  out_temp <- nleqslv(x= c(0.5,1,1), fn = score_func,control=list(ftol=1e-14,xtol=1e-14))
  pr_X_Y0Z0 <- out_temp$x[1]
  rr_X_Y_Z <- out_temp$x[2]
  rr_X_Z_Y <- out_temp$x[3]
  pr_X_Y1 <- pr_X_Y0 * rr_X_Y
  pr_X_Y1Z0 <- pr_X_Y0Z0 * rr_X_Y_Z
  pr_X_Y0Z1 <- pr_X_Y0Z0 * rr_X_Z_Y
  pr_X_Y1Z1 <- pr_X_Y0Z0 * rr_X_Z_Y * rr_X_Y_Z
  stopifnot(0<pr_X_Y0Z0 & pr_X_Y0Z0<1)
  stopifnot(isTRUE(all.equal(pr_X_Y0, pr_Z_Y0 * pr_X_Y0Z1 + (1-pr_Z_Y0) * pr_X_Y0Z0))) #beta0
  stopifnot(isTRUE(all.equal(pr_X_Y1, pr_Z_Y1 * pr_X_Y1Z1 + (1-pr_Z_Y1) * pr_X_Y1Z0)))
  pr_X_Z0 <- (1 - pr_Y_Z0) * pr_X_Y0Z0 + pr_Y_Z0 * pr_X_Y1Z0
  pr_X_Z1 <- (1 - pr_Y_Z1) * pr_X_Y0Z1 + pr_Y_Z1 * pr_X_Y1Z1
  stopifnot(isTRUE(all.equal(rr_X_Z, pr_X_Z1/pr_X_Z0)))
  pr_X <- pr_Z * pr_X_Z1 + (1-pr_Z) * pr_X_Z0
  stopifnot(isTRUE(all.equal(pr_X, pr_Y * pr_X_Y1 + (1-pr_Y) * pr_X_Y0)))
  
  pr_Z_X1 <- (pr_X_Z1 * pr_Z) / pr_X
  pr_Z_X0 <- ((1-pr_X_Z1) * pr_Z) / (1-pr_X)
  
  #using the structure of the DAG to encode this independence X indep A given Y,Z
  pr_X_Y0Z0A0 <- pr_X_Y0Z0A1 <- pr_X_Y0Z0
  pr_X_Y1Z0A0 <- pr_X_Y1Z0A1 <- pr_X_Y1Z0
  pr_X_Y0Z1A0 <- pr_X_Y0Z1A1 <- pr_X_Y0Z1
  pr_X_Y1Z1A0 <- pr_X_Y1Z1A1 <- pr_X_Y1Z1
  
  #a few last preliminary quantities
  pr_X_A1Z1 <- pr_X_Y1Z1A1 * pr_Y_A1Z1 + pr_X_Y0Z1A1 * (1-pr_Y_A1Z1)
  pr_X_A1Z0 <- pr_X_Y1Z0A1 * pr_Y_A1Z0 + pr_X_Y0Z0A1 * (1-pr_Y_A1Z0)
  pr_X_A0Z1 <- pr_X_Y1Z1A0 * pr_Y_A0Z1 + pr_X_Y0Z1A0 * (1-pr_Y_A0Z1)
  pr_X_A0Z0 <- pr_X_Y1Z0A0 * pr_Y_A0Z0 + pr_X_Y0Z0A0 * (1-pr_Y_A0Z0)
  pr_A_Z0 <- pr_A / (pr_Z * (rr_A_Z-1) + 1)
  pr_A_Z1 <- pr_A_Z0 * rr_A_Z
  
  #now, the complete marginal pmf of the observables
  pr_X1A1Z1 <- pr_X_A1Z1 * pr_A_Z1 * pr_Z
  pr_X0A1Z1 <- (1-pr_X_A1Z1) * pr_A_Z1 * pr_Z
  pr_X1A1Z0 <- pr_X_A1Z0 * pr_A_Z0 * (1-pr_Z)
  pr_X0A1Z0 <- (1-pr_X_A1Z0) * pr_A_Z0 * (1-pr_Z)
  pr_X1A0Z1 <- pr_X_A0Z1 * (1-pr_A_Z1) * pr_Z
  pr_X1A0Z0 <- pr_X_A0Z0 * (1-pr_A_Z0) * (1-pr_Z)
  pr_X0A0Z1 <- (1-pr_X_A0Z1) * (1-pr_A_Z1) * pr_Z
  pr_X0A0Z0 <- (1-pr_X_A0Z0) * (1-pr_A_Z0) * (1-pr_Z)
  stopifnot(isTRUE(all.equal(1,sum(pr_X1A1Z1, pr_X0A1Z1,
                                   pr_X1A1Z0, pr_X0A1Z0,
                                   pr_X1A0Z1, pr_X1A0Z0,
                                   pr_X0A0Z1, pr_X0A0Z0))))
  stopifnot(isTRUE(all.equal(pr_A, pr_X1A1Z1 + pr_X0A1Z1 + pr_X1A1Z0 + pr_X0A1Z0)))
  stopifnot(isTRUE(all.equal(pr_X, pr_X1A1Z1 + pr_X1A0Z1 + pr_X1A1Z0 + pr_X1A0Z0)))
  stopifnot(isTRUE(all.equal(pr_Z, pr_X1A1Z1 + pr_X1A0Z1 + pr_X0A1Z1 + pr_X0A0Z1)))
  
  #and the corresponding conditional probabilities viewing A as the outcome
  pr_A_X0Z0 <- pr_X0A1Z0 / ( (1-pr_X_Z0) * (1-pr_Z) )
  pr_A_X1Z0 <- pr_X1A1Z0 / ( (pr_X_Z0) * (1-pr_Z) )
  pr_A_X0Z1 <- pr_X0A1Z1 / ( (1-pr_X_Z1) * pr_Z )
  pr_A_X1Z1 <- pr_X1A1Z1 / ( (pr_X_Z1) * pr_Z )
  #these are the logistic regression coefficients under an interaction model
  beta0_int <- qlogis(pr_A_X0Z0)
  beta1_int <- qlogis(pr_A_X1Z0)-qlogis(pr_A_X0Z0)
  beta2_int <- qlogis(pr_A_X0Z1)-qlogis(pr_A_X0Z0)
  beta3_int <- (qlogis(pr_A_X1Z1)-qlogis(pr_A_X1Z0)) - 
    (qlogis(pr_A_X0Z1)-qlogis(pr_A_X0Z0))
  
  #finally, here is the objective function corresponding to the expected log likelihood
  #of the main effects logistic regression model for A with covariaets X and Z, taken over
  #the pmf of the true data generating mechanism
  #these get the 'true' coefficients from a logistic regression model adjusting for Z
  
  obj_func <- function(x,weighted="none",adj = TRUE){
    beta0 <- x[1]
    beta1 <- x[2]
    if(adj) beta2 <- x[3] else beta2 <- 0
    # beta3 <- x[4]
    
    #now, we're taking expectation over X,A,Z based on "true" DGM
    #of likelihood corresponding to a main effects logistic regression of A on X and Z
    ll_func <- function(X,A,Z){
      p <- plogis(q = beta0 + beta1 * X + beta2 * Z)
      # p <- plogis(q = beta0 + beta1 * X + beta2 * Z + beta3 * Z * X)
      A * log(p) + (1-A) * log1p(-p)
    }
    if(weighted=="none"){
      w_X1A1Z1 <- w_X1A1Z0 <- w_X1A0Z1 <- w_X1A0Z0 <-
        w_X0A1Z1 <- w_X0A1Z0 <- w_X0A0Z1 <- w_X0A0Z0 <- 1
    } else if(weighted=="ipw"){
      w_X1A1Z1 <- w_X1A0Z1 <- 1/pr_X_Z1
      w_X0A1Z1 <- w_X0A0Z1 <- 1/(1-pr_X_Z1)
      w_X1A1Z0 <- w_X1A0Z0 <- 1/pr_X_Z0
      w_X0A1Z0 <- w_X0A0Z0 <- 1/(1-pr_X_Z0)
    } else if(weighted=="ipw_att"){
      w_X1A1Z1 <- w_X1A0Z1 <- w_X1A1Z0 <- w_X1A0Z0 <- 1
      w_X0A1Z1 <- w_X0A0Z1 <- pr_X_Z1/(1-pr_X_Z1)
      w_X0A1Z0 <- w_X0A0Z0 <- pr_X_Z0/(1-pr_X_Z0)
    } else if(weighted=="bal"){
      w_X1A1Z1 <- w_X0A1Z1 <- 1/pr_A_Z1
      w_X1A1Z0 <- w_X0A1Z0 <- 1/pr_A_Z0
      w_X1A0Z1 <- w_X0A0Z1 <- 1/(1-pr_A_Z1)
      w_X1A0Z0 <- w_X0A0Z0 <- 1/(1-pr_A_Z0)
    } else if(weighted=="bal_att"){
      w_X1A1Z1 <- w_X0A1Z1 <- w_X1A1Z0 <- w_X0A1Z0 <- 1
      w_X1A0Z1 <- w_X0A0Z1 <- pr_A_Z1/(1-pr_A_Z1)
      w_X1A0Z0 <- w_X0A0Z0 <- pr_A_Z0/(1-pr_A_Z0)
    }
    ell_out <- 
      w_X1A1Z1 * pr_X1A1Z1 * ll_func(X=1,A=1,Z=1) + 
      w_X0A1Z1 * pr_X0A1Z1 * ll_func(X=0,A=1,Z=1) + 
      w_X1A1Z0 * pr_X1A1Z0 * ll_func(X=1,A=1,Z=0) + 
      w_X0A1Z0 * pr_X0A1Z0 * ll_func(X=0,A=1,Z=0) + 
      w_X1A0Z1 * pr_X1A0Z1 * ll_func(X=1,A=0,Z=1) + 
      w_X1A0Z0 * pr_X1A0Z0 * ll_func(X=1,A=0,Z=0) + 
      w_X0A0Z1 * pr_X0A0Z1 * ll_func(X=0,A=0,Z=1) + 
      w_X0A0Z0 * pr_X0A0Z0 * ll_func(X=0,A=0,Z=0)
    
    -ell_out
  }
  
  #these are 'numerically' estimated outputs after main effect logistic regression adjustment
  adj_temp <- TRUE
  betas_out <- optim(par=numeric(2 + adj_temp), #par=c(0,0,0,0), 
                     fn=obj_func, method = "BFGS", weighted="none", adj=adj_temp,
                     control=list(maxit=1000,reltol=1e-12))$par
  OR_adj <- exp(betas_out[2])
  betas_out <- optim(par=numeric(2 + adj_temp), #par=c(0,0,0,0), 
                     fn=obj_func, method = "BFGS", weighted="ipw", adj=adj_temp,
                     control=list(maxit=1000,reltol=1e-12))$par
  OR_IPW_adj <- exp(betas_out[2])
  betas_out <- optim(par=numeric(2 + adj_temp), #par=c(0,0,0,0), 
                     fn=obj_func, method = "BFGS", weighted="ipw_att", adj=adj_temp,
                     control=list(maxit=1000,reltol=1e-12))$par
  OR_IPWATT_adj <- exp(betas_out[2])
  betas_out <- optim(par=numeric(2 + adj_temp), #par=c(0,0,0,0), 
                     fn=obj_func, method = "BFGS", weighted="bal", adj=adj_temp,
                     control=list(maxit=1000,reltol=1e-12))$par
  OR_BAL_adj <- exp(betas_out[2])
  betas_out <- optim(par=numeric(2 + adj_temp), #par=c(0,0,0,0), 
                     fn=obj_func, method = "BFGS", weighted="bal_att", adj=adj_temp,
                     control=list(maxit=1000,reltol=1e-12))$par
  OR_BALATT_adj <- exp(betas_out[2])
  
  #marginal unadjusted OR
  pr_X_A1 <- pr_X_A1Z1 * pr_Z_A1 + pr_X_A1Z0 * (1-pr_Z_A1)
  pr_X_A0 <- pr_X_A0Z1 * pr_Z_A0 + pr_X_A0Z0 * (1-pr_Z_A0)
  stopifnot(isTRUE(all.equal(pr_X, pr_A * pr_X_A1 + (1-pr_A) * pr_X_A0)))
  OR_unadj <- pr_X_A1/(1-pr_X_A1) / ( pr_X_A0 / (1-pr_X_A0) )
  
  #marginal adjusted OR (balancing weights overall)
  pr_X_A1_BAL <- pr_X_A1Z1 * pr_Z + pr_X_A1Z0 * (1-pr_Z)
  pr_X_A0_BAL <- pr_X_A0Z1 * pr_Z + pr_X_A0Z0 * (1-pr_Z)
  OR_BAL <- pr_X_A1_BAL/(1-pr_X_A1_BAL) / ( pr_X_A0_BAL / (1-pr_X_A0_BAL) )
  #marginal adjusted OR (balancing weights to the infected)
  pr_X_A1_BALATT <- pr_X_A1Z1 * pr_Z_A1 + pr_X_A1Z0 * (1-pr_Z_A1)
  pr_X_A0_BALATT <- pr_X_A0Z1 * pr_Z_A1 + pr_X_A0Z0 * (1-pr_Z_A1)
  OR_BALATT <- pr_X_A1_BALATT/(1-pr_X_A1_BALATT) / ( pr_X_A0_BALATT / (1-pr_X_A0_BALATT) )
  
  #marginal unadjusted OR (reversed)
  pr_A_X1 <- pr_A_X1Z1 * pr_Z_X1 + pr_A_X1Z0 * (1-pr_Z_X1)
  pr_A_X0 <- pr_A_X0Z1 * pr_Z_X0 + pr_A_X0Z0 * (1-pr_Z_X0)
  OR_unadj3 <- pr_A_X1/(1-pr_A_X1) / ( pr_A_X0 / (1-pr_A_X0) )
  #marginal adjusted OR (reversed, IPW)
  pr_A_X1_IPW <- pr_A_X1Z1 * pr_Z + pr_A_X1Z0 * (1-pr_Z)
  pr_A_X0_IPW <- pr_A_X0Z1 * pr_Z + pr_A_X0Z0 * (1-pr_Z)
  OR_IPW <- pr_A_X1_IPW/(1-pr_A_X1_IPW) / ( pr_A_X0_IPW / (1-pr_A_X0_IPW) )
  #marginal adjusted OR (reversed, IPW to those with the symptom)
  pr_A_X1_IPWATT <- pr_A_X1Z1 * pr_Z_X1 + pr_A_X1Z0 * (1-pr_Z_X1)
  pr_A_X0_IPWATT <- pr_A_X0Z1 * pr_Z_X1 + pr_A_X0Z0 * (1-pr_Z_X1)
  OR_IPWATT <- pr_A_X1_IPWATT/(1-pr_A_X1_IPWATT) / ( pr_A_X0_IPWATT / (1-pr_A_X0_IPWATT) )
  
  pr_Y_A1 #pi
  pr_A #alpha
  #finally, compute the "theoretical" OR under no confounding
  pr_X_A1_nc <- ( rr_X_Y + (1-pr_Y_A1)/pr_Y_A1 ) / 
    ( rr_X_Y + (1-pr_A*pr_Y_A1)/(pr_A*pr_Y_A1) )
  pr_X_A0_nc <- ( (1-pr_X_Y0*rr_X_Y)*pr_Y_A1 + (1-pr_X_Y0)*(1-pr_Y_A1) ) / 
    ( (1-pr_X_Y0*rr_X_Y)*pr_Y_A1 + (1-pr_X_Y0)*(1/pr_A-pr_Y_A1) )
  OR_nc <- pr_X_A1_nc/(1-pr_X_A1_nc) / ( pr_X_A0_nc / (1-pr_X_A0_nc) )
  
  c(pr_A=pr_A, #alpha
    pr_Y_A1=pr_Y_A1, #pi
    pr_X_Y0=pr_X_Y0, #beta0
    rr_X_Y=rr_X_Y, #beta1
    pr_Z=pr_Z,
    pr_A_Z1=pr_A_Z1,pr_A_Z0=pr_A_Z0,
    pr_Z_A1=pr_Z_A1,pr_Z_A0=pr_Z_A0,
    rr_A_Z=rr_A_Z,
    rr_X_Z=rr_X_Z,
    pr_X_Y0Z0=pr_X_Y0Z0, rr_X_Y_Z=rr_X_Y_Z, rr_X_Z_Y=rr_X_Z_Y,
    pr_X_Z0=pr_X_Z0,
    pr_Y_A1Z0=pr_Y_A1Z0,pr_Y_A1Z1=pr_Y_A1Z1,
    pr_X_A1=pr_X_A1, pr_X_A0=pr_X_A0, 
    pr_X_A1_nc=pr_X_A1_nc,  pr_X_A0_nc=pr_X_A0_nc,  
    pr_A_X1=pr_A_X1, pr_A_X0=pr_A_X0, 
    OR_nc=OR_nc, OR_unadj=OR_unadj, 
    OR_unadj3=OR_unadj3, OR_adj=OR_adj,
    pr_A_X1_IPW=pr_A_X1_IPW, pr_A_X0_IPW=pr_A_X0_IPW,
    pr_A_X1_IPWATT=pr_A_X1_IPWATT, pr_A_X0_IPWATT=pr_A_X0_IPWATT,
    OR_IPW=OR_IPW,
    OR_IPWATT=OR_IPWATT,
    OR_IPW_adj=OR_IPW_adj,OR_IPWATT_adj=OR_IPWATT_adj,
    pr_X_A1_BAL=pr_X_A1_BAL,pr_X_A0_BAL=pr_X_A0_BAL,
    pr_X_A1_BALATT=pr_X_A1_BALATT,pr_X_A0_BALATT=pr_X_A0_BALATT,
    OR_BAL=OR_BAL,
    OR_BALATT=OR_BALATT,
    OR_BAL_adj=OR_BAL_adj,OR_BALATT_adj=OR_BALATT_adj)
}



cairo_pdf(filename = "OR_RR_conf_plots.pdf",
          width = 7,height=7)
par(mfrow=c(2,2),oma=c(4,2,0,0)) #default oma=c(0,0,0,0)
setting_mat <- NULL
for(i in 1:4){
  beta0_temp <- 1/5
  pi_temp <- 1/4
  alpha_temp <- 0.7
  pr_Z_temp <- 0.55
  rr_A_Z_temp <- c(4/3,3/4,2,0.5)
  rr_X_Z_temp <- c(4/3,4/3,2,2)

  beta1_max <- 1/beta0_temp
  beta1_min <- 0
  beta1 <- c(seq(from=beta1_min,to=1,length.out = 25)[-1],
             seq(from=1,to=beta1_max,length.out = 25)[-1])
  phi <- (1-pi_temp)/pi_temp * (1-beta0_temp)/beta0_temp
  out_mat <- t(sapply(beta1, function(x) OR_func(pr_A=alpha_temp,pr_Y_A1 = pi_temp,
                                                 pr_X_Y0 = beta0_temp, rr_X_Y = x,
                                                 pr_Z=pr_Z_temp, rr_A_Z = rr_A_Z_temp[i],
                                                 rr_X_Z = rr_X_Z_temp[i])))
  
  # palette <- RColorBrewer::brewer.pal(n=6,name = "Paired")
  palette <- RColorBrewer::brewer.pal(n=2,name = "Paired")
  matplot(x=beta1,y=out_mat[,c("OR_nc",
                               "OR_unadj",#"OR_adj",
                               "OR_BALATT",#"OR_BALATT_adj",
                               NULL)],type="l",
          col=c("grey80",palette),lwd=2,lty=1,
          xlim = c(0.15,1/0.15), ylim = c(0.15,1/0.15), log="xy",
          xlab="",ylab="",cex = 0.7,
          main=bquote(.(LETTERS[i]):~~~RR[paste(A,",",Z)]==.(round(rr_A_Z_temp[i],3))~~~~~~~~~RR[paste(X,",",Z)]==.(round(rr_X_Z_temp[i],3))))
  abline(a=0,b=1,lty=3,lwd=2)
  # abline(v=phi,lty=3,lwd=1,col="black")
  abline(v=1,lty=3,lwd=1,col="grey85")
  abline(h=1,lty=3,lwd=1,col="grey85")
  if(i==3){
    legend(x="bottomright", 
           legend=c(expression(OR[paste(A,",",X)]:~No~Confounding),
                    expression(OR[paste(A,",",X)]:~Unweighted),
                    expression(OR[paste(A,",",X)]^(w):~Bal.~Weights)),
           fill=c("grey80",palette),cex=0.7)
  }
  setting_mat <- rbind(setting_mat,out_mat[1,])
}
mtext(expression(True~Risk~Ratio~of~X~by~PASC~(beta[1])),
      side=1,line=0,outer=TRUE,cex=1.3)
mtext(expression(Analysis~Odds~Ratio~of~A~by~X),
      side=2,line=0,outer=TRUE,cex=1.3,las=0)
mtext(bquote(pi==.(round(pi_temp,3))~~~beta[0]==.(round(beta0_temp,3))~~~alpha==.(round(alpha_temp,3))~~~Pr(Z==1)==.(round(pr_Z_temp,3))),
      side=1,line=1,outer=TRUE,cex=0.9)
par(mfrow=c(1,1),oma=c(0,0,0,0))

dev.off()




