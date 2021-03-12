library(PerformanceAnalytics)
library(mcmc)
library(coda)
data(logit)

# Importing data
x<-read.table("/Users/manoj/desktop/7089CEM - Introduction to Statistical methods/Course Work/x.csv",header=FALSE,sep=",",dec=".")
y<-read.table("/Users/manoj/desktop/7089CEM - Introduction to Statistical methods/Course Work/y.csv",header=FALSE,sep=",",dec=".")
time<-read.table("/Users/manoj/desktop/7089CEM - Introduction to Statistical methods/Course Work/time.csv",header=FALSE,sep=",",dec=".")
data<-data.frame(time,y,x)
names(data)<-c("time","y","x1","x2","x3","x4")
head(data)

#                        Task 1: Preliminary data analysis

# Time series plots (input and output EEG signals)

# Output Signal
plot(data$time,data$y, type="l", col="red")

# Input Signal
par(mfrow=c(2,2))
plot(data$time,data$x1, type="l", col="red",main="EEG signal x1")
plot(data$time,data$x2, type="l", col="red",main="EEG signal x2")
plot(data$time,data$x3, type="l", col="red",main="EEG signal x3")
plot(data$time,data$x4, type="l", col="red",main="EEG signal x4")

# Distribution for each EEG signal

hist(data$y,main="Distribution of output signal y")

par(mfrow=c(2,2))
hist(data$x1,main="Distribution of signal x1")
hist(data$x2,main="Distribution of signal x2")
hist(data$x3,main="Distribution of signal x3")
hist(data$x4,main="Distribution of signal x4")

# Correlation and scatter plots (between different input EEG signals and the output EEG)

chart.Correlation(data[,-1], histogram=TRUE, pch=19)

#          Task 2: Regression- modelling the relationship between EEG signals

library(PerformanceAnalytics)
library(mcmc)
library(coda)
data(logit)

x<-read.table("/Users/manoj/desktop/7089CEM - Introduction to Statistical methods/Course Work/x.csv",header=FALSE,sep=",",dec=".")
y<-read.table("/Users/manoj/desktop/7089CEM - Introduction to Statistical methods/Course Work/y.csv",header=FALSE,sep=",",dec=".")
time<-read.table("/Users/manoj/desktop/7089CEM - Introduction to Statistical methods/Course Work/time.csv",header=FALSE,sep=",",dec=".")
data<-data.frame(time,y,x)
names(data)<-c("time","y","x1","x2","x3","x4")

model1<-data.frame(z1=data$x4,z2=data$x1^2,z3=data$x1^3,z4=data$x3^4,z0=1)
head(model1)

model2<-data.frame(z1=data$x3^3,z2=data$x3^4,z0=1)
head(model2)

model3<-data.frame(z1=data$x2,z2=data$x1^3,z3=data$x3^4,z0=1)
head(model3)

model4<-data.frame(z1=data$x4,z2=data$x1^3,z3=data$x3^4,z0=1) 
head(model4)

model5<-data.frame(z1=data$x4,z2=data$x1^2,z3=data$x1^3,z4=data$x3^4,z5=data$x1^4,z0=1) 
head(model5)

# Task 2.1

ols<-function(X,y)
     {
       X<-as.matrix(X)
       Y<-as.matrix(y)
       X.X<-solve(t(X) %*% X)       
       b<-b<-X.X %*% (t(X) %*% Y)
       colnames(b)<- "OLS_parm" 
       return(b)}
X1<-as.matrix(model1)
pram_m1<-solve(t(X1) %*% X1) %*% (t(X1) %*% data$y)
# Model-1 parameters 
pram_m1

ols<-function(X,y)
     {
       X<-as.matrix(X)
       Y<-as.matrix(y)
       X.X<-solve(t(X) %*% X)       
       b<-b<-X.X %*% (t(X) %*% Y)
       colnames(b)<- "OLS_parm" 
       return(b)}
X2<-as.matrix(model2)
pram_m2<-solve(t(X2) %*% X2) %*% (t(X2) %*% data$y)
# Model-2 parameters 
pram_m2

ols<-function(X,y)
     {
       X<-as.matrix(X)
       Y<-as.matrix(y)
       X.X<-solve(t(X) %*% X)       
       b<-b<-X.X %*% (t(X) %*% Y)
       colnames(b)<- "OLS_parm" 
       return(b)}
X3<-as.matrix(model3)
pram_m3<-solve(t(X3) %*% X3) %*% (t(X3) %*% data$y)
# Model-3 parameters 
pram_m3

ols<-function(X,y)
     {
       X<-as.matrix(X)
       Y<-as.matrix(y)
       X.X<-solve(t(X) %*% X)       
       b<-b<-X.X %*% (t(X) %*% Y)
       colnames(b)<- "OLS_parm" 
       return(b)}
X4<-as.matrix(model4)
pram_m4<-solve(t(X4) %*% X4) %*% (t(X4) %*% data$y)
# Model-4 parameters 
pram_m4

ols<-function(X,y)
     {
       X<-as.matrix(X)
       Y<-as.matrix(y)
       X.X<-solve(t(X) %*% X)       
       b<-b<-X.X %*% (t(X) %*% Y)
       colnames(b)<- "OLS_parm" 
       return(b)}
X5<-as.matrix(model5)
pram_m5<-solve(t(X5) %*% X5) %*% (t(X5) %*% data$y)
# Model-5 parameters 
pram_m5

# Task 2.2

rss_m1<-sum((data$y - X1 %*% pram_m1)^2)
rss_m2<-sum((data$y - X2 %*% pram_m2)^2)
rss_m3<-sum((data$y - X3 %*% pram_m3)^2)
rss_m4<-sum((data$y - X4 %*% pram_m4)^2)
rss_m5<-sum((data$y - X5 %*% pram_m5)^2)
rss<-data.frame(Model=paste("Model",1:5, sep="-"),
                RSS=c(rss_m1,rss_m2,rss_m3,rss_m4,rss_m5))
rss

# Task 2.3

n<- dim(data)[1]
A<- -(n/2)*log(2*pi)
log_lik1<- A - (n/2)*log(rss_m1/(n-1)) -(1/2*rss_m1/(n-1))*rss_m1
log_lik2<- A - (n/2)*log(rss_m2/(n-1)) -(1/2*rss_m2/(n-1))*rss_m2
log_lik3<- A - (n/2)*log(rss_m3/(n-1)) -(1/2*rss_m3/(n-1))*rss_m3
log_lik4<- A - (n/2)*log(rss_m4/(n-1)) -(1/2*rss_m4/(n-1))*rss_m4
log_lik5<- A - (n/2)*log(rss_m5/(n-1)) -(1/2*rss_m5/(n-1))*rss_m5
likelihood<-data.frame(Model=paste("Model",1:5, sep="-"),
            Log_Likelihood=c(log_lik1,log_lik2,log_lik3,log_lik4,log_lik5))
likelihood

# Task 2.4 

aic_m1<- 2*length(pram_m1) - 2*log_lik1
aic_m2<- 2*length(pram_m2) - 2*log_lik2
aic_m3<- 2*length(pram_m3) - 2*log_lik3
aic_m4<- 2*length(pram_m4) - 2*log_lik4
aic_m5<- 2*length(pram_m5) - 2*log_lik5

bic_m1<- length(pram_m1)*log(n) - 2*log_lik1
bic_m2<- length(pram_m2)*log(n) - 2*log_lik2
bic_m3<- length(pram_m3)*log(n) - 2*log_lik3
bic_m4<- length(pram_m4)*log(n) - 2*log_lik4
bic_m5<- length(pram_m5)*log(n) - 2*log_lik5

aic_bic<-data.frame(Model=paste("Model",1:5, sep="-"),
            AIC=c(aic_m1,aic_m2,aic_m3,aic_m4,aic_m5),
            BIC=c(bic_m1,bic_m2,bic_m3,bic_m4,bic_m5))
aic_bic

# Task 2.5

e1<- data$y - X1 %*% pram_m1
e2<- data$y - X2 %*% pram_m2
e3<- data$y - X3 %*% pram_m3
e4<- data$y - X4 %*% pram_m4
e5<- data$y - X5 %*% pram_m5

# Error distribution
par(mfrow=c(3,2))
hist(e1,main="Distribution of residual in model-1")
hist(e2,main="Distribution of residual in model-2")
hist(e3,main="Distribution of residual in model-3")
hist(e4,main="Distribution of residual in model-4")
hist(e5,main="Distribution of residual in model-5")

# qq-plot for all model residual
par(mfrow=c(3,2))
qqnorm(e1, pch = 1, frame = FALSE,main="QQ-plot for Model-1")
qqline(e1, col = "steelblue", lwd = 2)

qqnorm(e2, pch = 1, frame = FALSE)
qqline(e2, col = "steelblue", lwd = 2)

qqnorm(e3, pch = 1, frame = FALSE)
qqline(e3, col = "steelblue", lwd = 2)

qqnorm(e4, pch = 1, frame = FALSE)
qqline(e4, col = "steelblue", lwd = 2)

qqnorm(e5, pch = 1, frame = FALSE)
qqline(e5, col = "steelblue", lwd = 2)


# Task 2.6

aic_m1<- 2*length(pram_m1) - 2*log_lik1
aic_m2<- 2*length(pram_m2) - 2*log_lik2
aic_m3<- 2*length(pram_m3) - 2*log_lik3
aic_m4<- 2*length(pram_m4) - 2*log_lik4
aic_m5<- 2*length(pram_m5) - 2*log_lik5

bic_m1<- length(pram_m1)*log(n) - 2*log_lik1
bic_m2<- length(pram_m2)*log(n) - 2*log_lik2
bic_m3<- length(pram_m3)*log(n) - 2*log_lik3
bic_m4<- length(pram_m4)*log(n) - 2*log_lik4
bic_m5<- length(pram_m5)*log(n) - 2*log_lik5

aic_bic<-data.frame(Model=paste("Model",1:5, sep="-"),
            AIC=c(aic_m1,aic_m2,aic_m3,aic_m4,aic_m5),
            BIC=c(bic_m1,bic_m2,bic_m3,bic_m4,bic_m5))
aic_bic

#The smallest AIC and BIC values was found at `Model-3` and the residual of this model are more normal comparative to others model. 

# Task 2.7

data_m3<-data.frame(t=data$t,model3, y=data$y) 
tr_sample <-sort(sample(n, 0.7*n))
tr_data<-data_m3[tr_sample ,]
ts_data<-data_m3[-tr_sample ,]

train_m3<- lm(y ~ z1+z2+z3 , data=tr_data)
test_fit <- predict(train_m3, newdata=ts_data,interval="confidence",level = 0.95, se.fit=TRUE)
res_df<- data.frame(ts_data,test_fit$fit,se=test_fit$se.fit)
head(res_df)

# Plot for prediction
plot(res_df$t,res_df$fit, type="l",col="blue")
arrows(x0=res_df$t, y0=res_df$fit-res_df$se, x1=res_df$t, y1=res_df$fit+res_df$se, code=3, angle=90, length=0.1)
points(res_df$t,res_df$y, pch=16)
lines(res_df$t,res_df$upr, type="l", col="red",lty=3)
lines(res_df$t,res_df$lwr, type="l", col="red",lty=3)

#                                       Task 3

two_pram<-sort(abs(pram_m3),decreasing =TRUE)[1:2]
two_pram  

# Choose uniform prior for theta_bias and theta_1
prior<- function () {
  theta_bias<-runif(1, 0, 1)
  theta_1<-runif(1, -0.5, 0.5)
  c(theta_bias,theta_1)
}


ABC_acceptance <- function(par,y,X, theta){

         parm_sel<-order(abs(theta),decreasing =TRUE)[1:2]
         theta[parm_sel]<-par
         y_sim<- as.matrix(X) %*% as.matrix(theta)
         diffmean <- abs(mean(y) - mean(y_sim))
      if(diffmean < 0.05) return(TRUE) else return(FALSE)
}


run_MCMC_ABC <- function(y,X, theta,iterations){
    chain = array(dim = c(iterations+1,2))
        chain[1,] = prior() 
    for (i in 1:iterations){     
        # proposalfunction
        proposal = prior() 
        
        if(ABC_acceptance(proposal,y,X,theta)){
            chain[i+1,] = proposal
        }else{
            chain[i+1,] = chain[i,]
        }
    }
    return(mcmc(chain))
}

posterior <- run_MCMC_ABC(data$y,model3, pram_m3,10000)
head(posterior)

summary(posterior)

plot(posterior)