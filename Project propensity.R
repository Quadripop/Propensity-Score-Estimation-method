list.of.packages = c("Matching","rpart","randomForest","gbm","twang","ipred","neuralnet",
                     "nnet","e1071","klaR","xtable","flexmix","AUC","Hmisc","Kendall","lattice")




#Importing libraries
library(tidyverse)
library(ggplot2)
library(rpart)
library(gbm)
library(ipred)
library(randomForest)
library(Matching)
library(twang)
library(caret)
library(e1071)
#Utility function: Sample correlation

sample.corr<- function(x,rho){
  y<-rho*((x-mean(x))/sqrt(var(x))) + sqrt(1-rho^2)*rnorm(length(x))
  return (y)
}
sim_data_df<-data.frame()

# function to generate the covariates
covfun<-function(size,scenarioT){
  w1<-rnorm(size,mean=0,sd=1)
  w2<-rnorm(size,mean=0,sd=1)
  w3<-rnorm(size,mean=0,sd=1)
  w4<-rnorm(size,mean=0,sd=1)
  w5<-sample.corr(w1,0.2)
  w6<-sample.corr(w2,0.9)
  w7<-rnorm(size,mean=0,sd=1)
  w8<-sample.corr(w3,0.2)
  w9<-sample.corr(w4,0.9)
  w10<-rnorm(size,mean=0,sd=1)
  
  # dichotomize some covariates to create binary variables
  w1<-ifelse(w1>mean(w1),1,0)
  w3<-ifelse(w3>mean(w3),1,0)
  w5<-ifelse(w5>mean(w5),1,0)
  w6<-ifelse(w6>mean(w6),1,0)
  w8<-ifelse(w8>mean(w8),1,0)
  w9<-ifelse(w9>mean(w9),1,0)
  
  #scenarios for data generating model according to Setoguchi et al
  
  # A: additivity and linearity
  # B: mild non-linearity
  # C: moderate non-linearity
  # D: mild non-additivity
  # E: mild non-additivity and non-linearity
  # F: moderate non-additivity
  # G: moderate non-additivity and non-linearity
  
  #Coefficients for treatment equation according to Setoguchi et al
  b0<-0
  b1<-0.8
  b2<--0.25
  b3<-0.6
  b4<--0.4
  b5<--0.8
  b6<--0.5
  b7<-0.7
  
  # Coefficients for outcome equation
  a0<--3.85
  a1<-0.3
  a2<--0.36
  a3<--0.73
  a4<--0.2
  a5<-0.71
  a6<--0.19
  a7<-0.26
  gamma1<--0.4
  
  #Scenarios for treatment assignment
  
  if (scenarioT=="A"){
    #model with additivity and linearity
    trueps<-(1+ exp(-(b0+b1*w1+b2*w2+b3*w3+b4*w4+b5*w5+b6*w6+b7*w7)))^-1
    }else 
      if(scenarioT=="B"){
      #model with mild non-linearity
    trueps<- (1+ exp(-(b0+b1*w1+b2*w2+b3*w3+b4*w4+b5*w5+b6*w6+b7*w7+b2*w2^2)))^-1
    }else 
      if(scenarioT=="C"){
      # model with moderate non-linearity
    trueps<-(1+ exp(-(b0+b1*w1+b2*w2+b3*w3+b4*w4+b5*w5+b6*w6+b7*w7+b2*w2^2
                      +b4*w4^2+b7*w7^2)))^-1
    }else 
      if(scenarioT=="D"){
        #model with mild non-additivity
      trueps<-(1+ exp(-(b0+b1*w1+b2*w2+b3*w3+b4*w4+b5*w5+b6*w6+b7*w7
          +b1*0.5*w1*w3+b2*0.7*w2*w4+b4*0.5*w4*w5+b5*0.5*w5*w6)))^-1
      }else
        if(scenarioT=="E"){
          #model with mild non-additivity and non-linearity
          trueps<-(1+ exp(-(b0+b1*w1+b2*w2+b3*w3+b4*w4+b5*w5+b6*w6+b7*w7
          + b2*w2^2+ b1*0.5*w1*w3+ b2*0.7*w2*w4+ b4*0.5*w4*w5+b5*0.5*w5*w6)))^-1
        }else
          if(scenarioT=="F"){
            #model with moderate non-additivity
            trueps<-(1+ exp(-(b0+b1*w1+b2*w2+b3*w3+b4*w4+b5*w5+b6*w6+b7*w7
              +b1*0.5*w1*w3+ b2*0.7*w2*w4+ b3*0.5*w3*w5+b4*0.7*w4*w6
              +b5*0.5*w5*w7+ b1*0.5*w1*w6 + b2*0.7*w2*w3 + b3*0.5*w3*w4 
              + b4*0.5*w4*w5 + b5*0.5*w5*w6  )))^-1
          }else{
            #scenario G
            #model with moderate non-additivity and non-linearity
            trueps<-(1+ exp(-(b0+b1*w1+b2*w2+b3*w3+b4*w4+b5*w5+b6*w6+b7*w7
                + b2*w2^2+ b4*w4^2 + b7*w7^2 + b1*0.5*w1*w3 + b2*0.7*w2*w4
                + b3*0.5*w3*w5 + b4*0.7*w4*w6 + b5*0.5*w5*w7 + b1*0.5*w1*w6
                + b2*0.7*w2*w3 + b3*0.5*w3*w4 + b4*0.5*w4*w5 + b5*0.5*w5*w6)))^-1
              }
  
  #binary treatment
  unif1<-runif(size,0,1)
  T<-ifelse(trueps>unif1,1,0)
  
  #Outcome model(continuous)
  #For scenario A-G
  Y<- a0 + a1*w1 + a2*w2 + a3*w3 + a4*w4 + a5*w8 + a6*w9 +a7*w10 + gamma1*T
  Y1<- a0 + a1*w1 + a2*w2 + a3*w3 + a4*w4 + a5*w8 + a6*w9 +a7*w10 + gamma1
  Y0<- a0 + a1*w1 + a2*w2 + a3*w3 + a4*w4 + a5*w8 + a6*w9 +a7*w10 
  indeff<-Y1-Y0
  
  #create simulation dataset
  
  sim_data<- as.data.frame(cbind(w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,T,Y,indeff))
 
  sim_data_df<<-sim_data
  
}


funsim<-function(x,psmethod,par="ATT"){
  #set.seed(123)  # for reproducibility
  #folds <- createFolds(x$T, k = 5, list = TRUE)
  #x$ps <- rep(NA, nrow(x)) 
  
  if(psmethod =="truelogit"){
    mod= glm(T~w1+w2+w3+w4, data=x, family=binomial)
    x$ps = mod$fitted
    
  } else if (psmethod =="logit"){
    #for (fold_idx in seq_along(folds)) {
      #train_idx <- setdiff(seq_len(nrow(x)), folds[[fold_idx]])
      #test_idx <- folds[[fold_idx]]
    
    #mod = glm(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,data=x[train_idx,],
           #  family=binomial)
    mod = glm(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,data=x,
                family=binomial)
    x$ps= mod$fitted
    #x$ps[test_idx] <- predict(mod, newdata = x[test_idx, ], type = "response")}
    #x$ps=(1+ exp(-(b0+b1*w1+b2*w2+b3*w3+b4*w4+b5*w5+b6*w6+b7*w7)))^-1
   
  } else if (psmethod =="tree"){
    mod = rpart(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,method="class",
                data=x,cp=0.01)
    x$ps = predict(mod,type="prob")[,2]
     #x$ps=(1+ exp(-(b0+b1*w1+b2*w2+b3*w3+b4*w4+b5*w5+b6*w6+b7*w7)))^-1

  } else if (psmethod =="randomforest"){
    #ctrl <- trainControl(method = "cv", number = 10)
    mod = randomForest(factor(T)~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,
                       data=x, ntree=10000,mtry=8,nodesize=10)
        #mod <- train(
             # factor(T) ~ ., data = x[, 1:11],
             # method = "rf",
              #trControl = ctrl,
              #tuneLength = 3  # try fewer or more trees
#)
    
    x$ps<- predict(mod,newdata=x, type = "prob")[,2]
  } else if (psmethod == "gbm") {
    ctrl <- trainControl(method = "cv", number = 5)
    #mod = gbm(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x , distribution = "bernoulli",
              #interaction.depth = 1, n.trees=1000)
    mod <- train(
      factor(T) ~ ., data = x[, 1:11],
      method = "gbm",
      trControl = ctrl,
      verbose = FALSE,
      tuneLength = 5
    )
    #x$ps = predict.gbm (mod, data=x , n.trees=1000, type="response")
    x$ps = predict(mod, newdata=x ,type="prob")[,2]
  } else if (psmethod =="gbmtwang"){
    mod = ps(T~w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x, n.trees=10000,
             interaction.depth=4, verbose= FALSE, shrinkage= 0.001,bag.fraction = 0.5)
    x$ps = as.vector(mod$ps[,1])
  } else if (psmethod == "bag"){
    mod = bagging(T~w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x)
    x$ps = predict(mod,newdata=x, type="prob")
  }
  
  
  #true ATT
  g<-mean(x$indeff[x$T==1])
  
  #true ATE
  if(par=="ATE"){g<-mean(x$indeff)}
  
  # Estimating ATT via propensity score weighting
  
  weights<- ifelse(x$T==1,1,x$ps/(1-x$ps))
  
  outcome_mod1 <- lm(Y ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data = x[x$T == 1, ])
  outcome_mod0 <- lm(Y ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data = x[x$T == 0, ])
  
  # Predict outcomes for all observations under both treatment and control
  mu1_hat <- predict(outcome_mod1, newdata = x)
  mu0_hat <- predict(outcome_mod0, newdata = x)
  
  dr_est <- mean(
    ((x$T - x$ps) / (x$ps * (1 - x$ps))) * (x$Y - ifelse(x$T == 1, mu1_hat, mu0_hat)) +
      mu1_hat - mu0_hat
  )
  
  #Estimating ATE via propensity score weighting
  
  if(par=="ATE"){
    
    x$ps <- pmin(pmax(x$ps, 0.01), 0.99) 
    weights<-ifelse(x$T==1,1/x$ps,1/(1-x$ps))
  
  p_weight<-hist(weights, breaks = 30, main = "Histogram of Weights")
  weight_5<-mean(weights>5)
  weight_10<-mean(weights>10)
  dr_absrbias <- abs((dr_est - g) / g) * 100
  dr_var <- (dr_est - g)^2}
  
  #if(par=="ATE"){# Truncate extreme propensity scores
    #x$ps <- pmin(pmax(x$ps, 0.01), 0.99)  # adjust  threshold if needed
    
    # Stabilized weights
    #p_treated <- mean(x$T == 1)
    #weights <- ifelse(x$T == 1, p_treated / x$ps, (1 - p_treated) / (1 - x$ps))
 # }
  g_hat<-weighted.mean(x$Y[x$T==1],weights=weights[x$T==1]) - 
    weighted.mean(x$Y[x$T==0],weights=weights[x$T==0])
  
  
  
  absrbias <- abs((g_hat-g)/g)*100
  varg_hat <- (g_hat-g)^2
  
  modw     <- lm( Y ~  T, data=x, weights=weights)
  hatgsew <- summary(modw)$coefficients[c("T"),c("Std. Error")]
  covw <- ifelse(g > g_hat-2*hatgsew  & g < g_hat + 2*hatgsew , 1, 0)
  
  ps_plot<-ggplot(x, aes(x = ps, fill = as.factor(T))) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
    labs(title = "Propensity Score Distribution", x = "Propensity Score", fill = "Treatment Status") +
    theme_minimal()
  ps_summary <- summary(x$ps)
  ps_extremes <- mean(x$ps < 0.1 | x$ps > 0.9)
  weight_summary <- summary(weights)
  
  # estimating ATT via matching with caliper
  
  rr = Match(Y=x$Y, Tr=x$T, X=x$ps, caliper=0.25, M=1, replace=TRUE,ties=FALSE)
  
  if(par=="ATE"){
    rr = Match(Y=x$Y, Tr=x$T, X=x$ps, caliper=0.25,M=1,estimand=par,replace=TRUE,ties=FALSE)
  }
  
  # estimating ATE via matching with caliper
  
  allmdata<-rbind(x[rr$index.treated,],x[rr$index.control,])
  
  hatgm          <-rr$est
  absrbiasm    <- abs((hatgm -g)/g)*100 
  varhatgm      <- (hatgm-g)^2 
  
  hatgsem          <-rr$se.standard
  covm <- ifelse(g > hatgm-2*hatgsem  & g < hatgm + 2*hatgsem , 1, 0)
  # = mean(mdata$Y[mdata$Tr==1])-mean(mdata$Y[mdata$Tr==0])
  
  # ~  size of matched dataset
  
  orig.nobs   <-rr$orig.nobs
  orig.tnobs  <-rr$orig.treated.nobs
  
  match.nobs   <-length(rr$mdata$Tr)
  match.tnobs  <-length(rr$mdata$Tr[rr$mdata$Tr==1])
  
  
  
  
  #balance before weighting for w1,w2,...w10
  
  bb = MatchBalance(T ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10
                    +w1*w2+w1*w3+w1*w4+w1*w5+w1*w6+w1*w7+w1*w8+w1*w9+w1*w10
                    +w2*w3+w2*w4+w2*w5+w2*w6+w2*w7+w2*w8+w2*w9+w2*w10
                    +w3*w4+w3*w5+w3*w6+w3*w7+w3*w8+w3*w9+w3*w10
                    +w4*w5+w4*w6+w4*w7+w4*w8+w4*w9+w4*w10
                    +w5*w6+w5*w7+w5*w8+w5*w9+w5*w10
                    +w6*w7+w6*w8+w6*w9+w6*w10
                    +w7*w8+w7*w9+w7*w10
                    +w8*w9+w8*w10
                    +w9*w10, data=x, ks = FALSE, nboots=0, print.level=0)
  
  
  
  # ASAM
  asb_b    <-vector();for(i in 1:10){asb_b[[i]] <- bb$BeforeMatching[[i]]$sdiff}
  
  # ASAM with interactions
  interasbb<-vector();for(i in 1:(length(bb$BeforeMatching))){interasbb[[i]] <- bb$BeforeMatching[[i]]$sdiff}
  
  #~~ performance metrics: balance AFTER weighting for w1 ...w10
  
  ba = MatchBalance(T ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10
                    +w1*w2+w1*w3+w1*w4+w1*w5+w1*w6+w1*w7+w1*w8+w1*w9+w1*w10
                    +w2*w3+w2*w4+w2*w5+w2*w6+w2*w7+w2*w8+w2*w9+w2*w10
                    +w3*w4+w3*w5+w3*w6+w3*w7+w3*w8+w3*w9+w3*w10
                    +w4*w5+w4*w6+w4*w7+w4*w8+w4*w9+w4*w10
                    +w5*w6+w5*w7+w5*w8+w5*w9+w5*w10
                    +w6*w7+w6*w8+w6*w9+w6*w10
                    +w7*w8+w7*w9+w7*w10
                    +w8*w9+w8*w10
                    +w9*w10, data=x,weights=weights, ks = FALSE, nboots=0, print.level=0)
  
  # ASAM
  asb_a    <-vector()
  for(i in 1:10){asb_a[[i]] <- ba$BeforeMatching[[i]]$sdiff}
  
  # ASAM with interactions 
  interasba<-vector()
  for(i in 1:(length(ba$BeforeMatching))){interasba[[i]] <- ba$BeforeMatching[[i]]$sdiff}
  
  #~~ performance metrics: balance AFTER matching for w1 ...w10
  
  bam = MatchBalance(T ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10
                     +w1*w2+w1*w3+w1*w4+w1*w5+w1*w6+w1*w7+w1*w8+w1*w9+w1*w10
                     +w2*w3+w2*w4+w2*w5+w2*w6+w2*w7+w2*w8+w2*w9+w2*w10
                     +w3*w4+w3*w5+w3*w6+w3*w7+w3*w8+w3*w9+w3*w10
                     +w4*w5+w4*w6+w4*w7+w4*w8+w4*w9+w4*w10
                     +w5*w6+w5*w7+w5*w8+w5*w9+w5*w10
                     +w6*w7+w6*w8+w6*w9+w6*w10
                     +w7*w8+w7*w9+w7*w10
                     +w8*w9+w8*w10
                     +w9*w10, data=allmdata, ks = FALSE, nboots=0, print.level=0) 
  
  # ASAM after matching
  asb_am    <-vector()
  for(i in 1:10){asb_am[[i]] <- bam$BeforeMatching[[i]]$sdiff}
  
  # ASAM with interactions after matching
  interasbam<-vector()
  for(i in 1:(length(bam$BeforeMatching))){interasbam[[i]] <- bam$BeforeMatching[[i]]$sdiff}
  
  
  output<- list("g_hat"=g_hat,
                "hatgm"=hatgm,
                "absrbias"=absrbias,
                "absrbiasm"=absrbiasm,
                "varg_hat"=varg_hat,
                "varhatgm"=varhatgm,
                "hatgsew"=hatgsew,
                "hatgsem"=hatgsem,
                "covw"=covw,
                "covm"=covm,
                #  mean (over covariates) of balance summaries : 
                # asam
                "masb_b"=mean(abs(asb_b)),
                "masb_aw"=mean(abs(asb_a)),
                "masb_am"=mean(abs(asb_am)),
                # asam with interactions
                "masbinter_b"=mean(abs(interasbb)),
                "masbinter_aw"=mean(abs(interasba)),
                "masbinter_am"=mean(abs(interasbam)),
                "ps_plot"=ps_plot,
                "ps_summary"=ps_summary,
                "ps_extremes"=ps_extremes,
                "weight_diagnostics"=weight_summary,
                "p_weight"=p_weight,
                "weight_5"=weight_5,
                "weight_10"=weight_10,
                "dr_est" = dr_est,
                "dr_absrbias" = dr_absrbias,
                "dr_var" = dr_var
                )
  return(output)
  
}







