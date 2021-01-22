
# Last edited on a Windows 10 machine, November 3, 2020.
# rbreslawski@smu.edu for questions.

# Load packages
library(rstan)
library(loo)
library(parallel)
library(matrixStats)

# Fix randomization
set.seed(342543)

# Import data
d <- read.csv("data_PHTfaunas.csv", header=TRUE, stringsAsFactors=FALSE)

# Remove any NA cases and create a column for all components
d <- d[!is.na(d$Paleol), ]

# Sort by extinction status then alphabetical genus
d <- d[with(d, order(Extinct, Genus)), ]
# Index values for each extinction status
d$E_ind <- sapply(1:nrow(d), function(x){
  length(which(d$Extinct[1:x]==d$Extinct[x]))
})
# Index all cases
d$I <- 1:nrow(d)
# Vector of extinction values: 1 = extant, 2 = extinct
d$E <- d$Extinct + 1

# Priors for Stan model parameters
phi_scale_prior <- 0.1
mu_alpha_prior <- 1
mu_beta_prior <- 1

# MCMC argument for Stan
adapt_delta <- 0.9

# Data for Stan models. Include values for priors.
Stan_data <- list(n=nrow(d), tn=nrow(d), Archaeol=d$Archaeol, 
                  Faunas=d$Paleol + d$Archaeol, 
                  E=d$E, Statuses=max(d$E),
                  I=d$I, ll_i=d$I, ll_n=nrow(d), 
                  phi_scale_prior=phi_scale_prior, 
                  mu_alpha_prior=mu_alpha_prior,
                  mu_beta_prior=mu_beta_prior)

# Stan model parameters
iter <- 1e4
warmup <- 5e3
chains <- 4

# Fit extinction model
model_Extinct <- stan("model_extinction.stan", data=Stan_data, 
                      iter=iter, warmup=warmup, chains=chains,
                      control=list(adapt_delta=adapt_delta))
summary_Extinct <- summary(model_Extinct)$summary

# Fit null model
model_Null <- stan("model_null.stan", data=Stan_data, 
                   iter=iter, warmup=warmup, chains=chains,
                   control=list(adapt_delta=adapt_delta))
summary_Null <- summary(model_Null)$summary

# Posterior samples
samples_Extinct <- extract(model_Extinct)
samples_Null <- extract(model_Null)

# Save data
save(d, samples_Extinct, samples_Null, model_Extinct, model_Null,
     summary_Extinct, summary_Null, Stan_data, file="PHT_Stan_data.RData")


####--- --- -- -  MODEL COMPARISON WITH LOO-CV  - -- --- ---###

# This section has three central parts. First, a function is 
# defined that fits Stan models to n datasets, where each dataset
# consists of a training sample and a hold out sample. After 
# fitting the models, the function extracts the log likelihoods
# for the hold out data in each fitted model and returns an 
# object summarizing these values for all held out observations
# across the n datasets.
# 
# The second part prepares datasets for the held out
# taxa. Finally, the function is applied to both models for
# both hold out scenarios, resulting in four objects that
# summarize log likelihood for held out observations. These
# results are saved to disk.

# PART 1.
# Function. Expected log pointwise predictive density (ELPD).
# This function takes a Stan file name string, a list of data 
# for the Stan models (one list object per data fold), iterations, 
# warmup, and chains arguments for Stan, an integer for number
# heldout per fold (h), and the number of folds (K). It returns a
# list with the pointwise predictive density across observations,
# the ELPD, and the ELPD standard error.
elpd_list <- function(m_name, s_data, iter, warmup, chains, h, K){
  
  # Load libraries if they were not loaded at the beginning
  # of the script.
  library(parallel)
  library(rstan)
  library(loo)
  library(matrixStats)
  
  cat(paste("\nCompiling", m_name, "...\n"))
  m <- stan_model(file=m_name)
  
  # Make cluster and export necessary data and packages
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, varlist=list("m", "iter", "warmup", 
                                 "chains", "adapt_delta"),
                envir=environment())
  clusterEvalQ(cl, library("rstan"))
  clusterEvalQ(cl, library("loo"))
  
  cat(paste("Fitting", length(s_data), "dataset folds...\n"))
  
  ll_list <- parLapply(cl, s_data, function(x){
    model_fit <- sampling(m, data=x, iter=iter, 
                          warmup=warmup, chains=chains,
                          control=list(adapt_delta=adapt_delta))
    return(extract_log_lik(model_fit))
  })
  
  stopCluster(cl)
  
  cat("Models fitted and log likelihoods extracted.\n")
  
  cols <- K*h
  rows <- nrow(ll_list[[1]])
  log_lik_m <- matrix(NA, ncol=cols, nrow=rows)
  for(i in 1:K){
    log_lik_m[ , ((i-1)*h + 1):(i*h)] <- ll_list[[i]]
  }
  
  cat("Log likelihood matrix complete.\n")
  
  # Function. Get log column means from
  # log likelihood matrix.
  logColMeansExp <- function(x) {
    samples <- nrow(x)
    return(colLogSumExps(x) - log(samples))
  }
  
  pointwise <-  matrix(logColMeansExp(log_lik_m), ncol=1)
  colnames(pointwise) <- "elpd"
  
  elpd <- sum(pointwise)
  se_elpd <-  sqrt(ncol(log_lik_m) * var(pointwise))
  
  elpd_list <- list(pointwise=pointwise, elpd=elpd, se_elpd=se_elpd)
  
  cat(paste(m_name, "ELPD complete.\n\n"))
  
  return(elpd_list)
}


# PART 2.
#### CROSS VALIDATION DATASETS CREATION: WITHOLD TAXA ####
#### 
### Create n datasets, where n is equal to the number of taxa
### in the complete sample. For each dataset, withhold all 
### observations for a taxon.
###
Data_list_WHT <- lapply(1:nrow(d), function(x){
  
  df <- rbind(d, d[x, ])
  df$Paleol[x] <- 0
  df$Archaeol[x] <- 0
  
  Stan_data <- list(n=nrow(df), tn=nrow(df) - 1, 
                    Archaeol=df$Archaeol, 
                    Faunas=df$Paleol + df$Archaeol, 
                    E=df$E, Statuses=max(df$E),
                    I=df$I, ll_n=1, 
                    ll_i=as.array(nrow(df)),
                    phi_scale_prior=phi_scale_prior, 
                    mu_alpha_prior=mu_alpha_prior,
                    mu_beta_prior=mu_beta_prior)
  return(Stan_data)
})

# PART 3.
# EXECUTE FUNCTION ACROSS BOTH DATASETS AND BOTH MODELS
# 
# Obtain ELPD for each model  with entire taxa withheld.
elpd_extinct_WHT <- elpd_list("model_extinction.stan", Data_list_WHT, 
                                iter, warmup, chains, 1, nrow(d))
elpd_null_WHT <- elpd_list("model_null.stan", Data_list_WHT, 
                             iter, warmup, chains, 1, nrow(d))
# Get model weights
model_weights_WHT <- stacking_weights(cbind(elpd_extinct_WHT$pointwise[, 1],
                                              elpd_null_WHT$pointwise[, 1]))

# Summarize ELPD for each withheld taxon.
ho_WHT <- d
ho_WHT$elpd_extinct <- elpd_extinct_WHT$pointwise
ho_WHT$elpd_null <- elpd_null_WHT$pointwise
ho_WHT$lik_extinct <- exp(ho_WHT$elpd_extinct)
ho_WHT$lik_null <- exp(ho_WHT$elpd_null)
ho_WHT$lik_ratio <- ho_WHT$lik_extinct/ho_WHT$lik_null
ho_WHT$log_lik_ratio <- log(ho_WHT$lik_ratio)

# Save model comparison data
save(elpd_extinct_WHT, elpd_null_WHT, model_weights_WHT, ho_WHT,
     file="PHT_model_comparison.RData")
