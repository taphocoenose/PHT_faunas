
# Last edited on a Windows 10 machine, November 10, 2020.
# rbreslawski@smu.edu for questions.

# Plotting packages
library(ggplot2)
library(patchwork)
library(HDInterval)

# Load data
load("PHT_Stan_data.RData")
load("PHT_model_comparison.RData")

# Logit and inverse logit functions
logit <- function(x) log(x/(1-x))
inv_logit <- function(x) exp(x)/(exp(x) + 1)

# Function: Estimate beta distribution pars 
# from samples and return density values
# over grid.
EstBeta <- function(x, grid){
  
  mu <- mean(x)
  var <- var(x)
  a <- ((1 - mu) / var - 1 / mu) * mu^2
  b <- a * (1 / mu - 1)
  
  df <- data.frame(x=grid,
                   dnsty=dbeta(grid, a, b))
  return(df)
}




####--- --- -- -  WRITE TABLES  - -- --- ---###

# Write tables used in the manuscript.

Table_1 <- data.frame(model=c(rep("Null", 2), rep("Extinction", 5)), 
                      par=c("mu", "phi", "mu_extant", "mu_extinct",
                            "phi_extant", "phi_extinct", "mu_diff"),
                      post_med=rep(NA, 7), post_2.5=rep(NA, 7), 
                      post_97.5=rep(NA, 7),
                      stringsAsFactors=FALSE)

Table_1$post_med[1] <- median(samples_Null$mu)
Table_1$post_med[2] <- median(samples_Null$phi)
Table_1$post_med[3:4] <- apply(samples_Extinct$mu, 2, median)
Table_1$post_med[5:6] <- apply(samples_Extinct$phi, 2, median)
Table_1$post_med[7] <- median(samples_Extinct$mu_diff_extinct_extant)

Table_1[1, 4:5] <- hdi(samples_Null$mu)
Table_1[2, 4:5] <- hdi(samples_Null$phi)
Table_1[3:4, 4:5] <- t(apply(samples_Extinct$mu, 2, hdi))
Table_1[5:6, 4:5] <- t(apply(samples_Extinct$phi, 2, hdi))
Table_1[7, 4:5] <- hdi(samples_Extinct$mu_diff_extinct_extant)

write.csv(Table_1, "Table_1.csv", row.names=FALSE)

Table_2 <- data.frame(Genus=ho_WHT$Genus,
                      Archaeol_occurrences=ho_WHT$Archaeol,
                      Paleol_occurrences=ho_WHT$Paleol,
                      Extinction_Status=ho_WHT$Extinct,
                      Null_postMed=rep(NA, nrow(ho_WHT)),
                      Null_post2.5=rep(NA, nrow(ho_WHT)),
                      Null_post97.5=rep(NA, nrow(ho_WHT)),
                      Null_ELPD=ho_WHT$elpd_null[, 1],
                      Extinct_postMed=rep(NA, nrow(ho_WHT)),
                      Extinct_post2.5=rep(NA, nrow(ho_WHT)),
                      Extinct_post97.5=rep(NA, nrow(ho_WHT)),
                      Extinct_ELPD=ho_WHT$elpd_null[, 1],
                      LogLikRatio=ho_WHT$log_lik_ratio[, 1],
                      stringsAsFactors=FALSE)

Table_2[, 5:7] <- t(sapply(1:nrow(d), function(x){
  c(median(samples_Null$Arch_prop[, x]), 
    hdi(samples_Null$Arch_prop[, x]))
}))

Table_2[, 9:11] <- t(sapply(1:nrow(d), function(x){
  c(median(samples_Extinct$Arch_prop[, x]), 
    hdi(samples_Extinct$Arch_prop[, x]))
}))

write.csv(Table_2, "Table_2.csv", row.names=FALSE)

Table_3 <- data.frame(Model=c("Extinct", "Null"),
                      ELPD=c(elpd_extinct_WHT$elpd,
                             elpd_null_WHT$elpd),
                      ELPDse=c(elpd_extinct_WHT$se_elpd,
                               elpd_null_WHT$se_elpd),
                      Weights=c(model_weights_WHT),
                      stringsAsFactors=FALSE)

write.csv(Table_3, "Table_3.csv", row.names=FALSE)


#####################################################################
#####################   FIGURE 1  ###################################
#####################################################################


# Create mu and phi values for which to plot beta densities
mus <- seq(0.1, 0.9, by=0.2)
phis <- c(1, 2, 5, 20, 100)
grid <- seq(0.005, 0.995, by=0.005)

for(i in 1:length(phis)){
  for(j in 1:length(mus)){
    
    y <- dbeta(grid, mus[j]*phis[i], (1-mus[j])*phis[i])
    ij_df <- data.frame(x=grid, y=y)
    
    b_plot <- ggplot(ij_df, aes(x=x, y=y))+
      scale_x_continuous(expand=c(0,0), 
                         breaks=c(0, 0.5, 1.0),
                         limits=c(0, 1))+
      scale_y_continuous(expand=c(0.02,0))+
      geom_line()+
      ylab(as.expression(bquote(italic(phi["average"])==.(phis[i]))))+
      xlab(as.expression(bquote(italic("p"["average"])==.(mus[j]))))
    
    if(i==length(phis) & j==1){
      
      b_plot <- b_plot +
        theme(panel.background=element_rect(fill="white", color="grey"),
              panel.grid=element_blank(), 
              axis.text.y=element_blank(),
              axis.text.x=element_text(angle=90, vjust=0.4),
              axis.ticks.y=element_blank(), 
              axis.ticks.x=element_line(color="grey"))
        
    }else if(j==1){
      
      b_plot <- b_plot +
        theme(panel.background=element_rect(fill="white", color="grey"),
              panel.grid=element_blank(), 
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank())

      
    }else if(i==length(phis)){
      
      b_plot <- b_plot +
        theme(panel.background=element_rect(fill="white", color="grey"),
              panel.grid=element_blank(), 
              axis.text.y=element_blank(),
              axis.text.x=element_text(angle=90, vjust=0.4),
              axis.title.y=element_blank(),
              axis.ticks.y=element_blank(), 
              axis.ticks.x=element_line(color="grey"))

      
    }else{
      
      b_plot <- b_plot +
        theme(panel.background=element_rect(fill="white", color="grey"),
              panel.grid=element_blank(), 
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              axis.title=element_blank())
      
    }
    
    if(i==1 & j==1){
      Figure_1 <- b_plot
    }else{
      Figure_1 <- Figure_1 + b_plot
    }
    
  }
  
}

Figure_1 <- Figure_1 + plot_layout(ncol=length(mus))

ggsave("Figure_1.jpeg", plot=Figure_1, device="jpeg", 
       dpi=1000, height=6, width=8, units="in")


#####################################################################
#####################   FIGURE 3  ###################################
#####################################################################

# Get vector positions for the order of modal
# posterior archaeological proportions.
prop_mu <- order(sapply(1:nrow(d),function(x){
  
  grid <- seq(0, 1, length.out=500)
  
  b <- EstBeta(samples_Null$Arch_prop[, x], grid)
  b <- b[!(b$dnsty%in%c(-Inf,Inf)), ]
  
  xpos <- b$x[which.max(b$dnsty)]
  return(xpos)
}))

# For each taxon, generating geometry to plot
# violin densities for posterior estimates of
# archaeological proportions. Two violin
# densities are created, one for the null model
# and one for the extinction model.
taxon_df <- lapply(1:length(prop_mu), function(y){
  
  x <- prop_mu[y]
  
  grid <- seq(0, 1, length.out=500)
  null_df <- EstBeta(samples_Null$Arch_prop[, x], grid)
  null_df$model <- rep("null", nrow(null_df))
  extinct_df <- EstBeta(samples_Extinct$Arch_prop[, x], grid)
  extinct_df$model <- rep("extinct", nrow(extinct_df))
  extinct_df$x <- extinct_df$x + 1.05
  
  plotdf <- rbind(null_df, extinct_df)
  plotdf$Genus <- rep(d$Genus[x], nrow(plotdf))
  plotdf$Extinct <- rep(d$Extinct[x], nrow(plotdf))
  
  plotdf <- plotdf[!(plotdf$dnsty%in%c(-Inf,Inf)), ]
  
  plotdf$y <- rep(y, nrow(plotdf))
  
  return(plotdf)
})
taxon_df <- do.call("rbind", taxon_df)
taxon_df$group <- paste0(taxon_df$Genus, taxon_df$model)

# Rescale y values for plotting so that densities are 
# proportional to eachother.
taxon_df$scale <- taxon_df$dnsty/max(taxon_df$dnsty)
taxon_df$ymin <- taxon_df$y - taxon_df$scale/2
taxon_df$ymax <- taxon_df$y + taxon_df$scale/2
# Remove very small densities for plotting
taxon_df <- taxon_df[which(taxon_df$scale > 1/200), ]

# posterior 95% HDIs and modal values for the average 
# genus in each model
modal_val <- function(x){
  de <- density(x, n=1e4)
  xpos <- which.max(de$y)
  return(de$x[xpos])
}
null_genus <- rep(NA, 3)
null_genus[1] <- modal_val(samples_Null$mu)
null_genus[2] <- hdi(samples_Null$mu)[[1]]
null_genus[3] <- hdi(samples_Null$mu)[[2]]
extinct_genus_extinct <- rep(NA, 3)
extinct_genus_extinct[1] <- modal_val(samples_Extinct$mu[ , 2])
extinct_genus_extinct[2] <- hdi(samples_Extinct$mu[ , 2])[[1]]
extinct_genus_extinct[3] <- hdi(samples_Extinct$mu[ , 2])[[2]]
extinct_genus_extant <- rep(NA, 3)
extinct_genus_extant[1] <- modal_val(samples_Extinct$mu[ , 1])
extinct_genus_extant[2] <- hdi(samples_Extinct$mu[ , 1])[[1]]
extinct_genus_extant[3] <- hdi(samples_Extinct$mu[ , 1])[[2]]


# Create plot by taxon
taxon_plot <- ggplot(taxon_df, aes(x=x, ymin=ymin, ymax=ymax))+
  annotate("rect", xmin=null_genus[2], xmax=null_genus[3], 
           ymin=0.5, ymax=nrow(d)+0.5, fill="purple", alpha=0.4)+
  annotate("segment", x=null_genus[1], xend=null_genus[1], size=1,
           y=0.5, yend=nrow(d)+0.5, color="white", linetype="dashed")+
  annotate("rect", xmin=extinct_genus_extinct[2]+1.05, 
           xmax=extinct_genus_extinct[3]+1.05, 
           ymin=0.5, ymax=nrow(d)+0.5, 
           fill="red", alpha=0.4)+
  annotate("segment", x=extinct_genus_extinct[1]+1.05, 
           xend=extinct_genus_extinct[1]+1.05, 
           size=1, y=0.5, yend=nrow(d)+0.5, color="white", linetype="dashed")+
  annotate("rect", xmin=extinct_genus_extant[2]+1.05, 
           xmax=extinct_genus_extant[3]+1.05, 
           ymin=0.5, ymax=nrow(d)+0.5, fill="blue", alpha=0.4)+
  annotate("segment", x=extinct_genus_extant[1]+1.05, 
           xend=extinct_genus_extant[1]+1.05, 
           size=1, y=0.5, yend=nrow(d)+0.5, color="white", linetype="dashed")+
  annotate("rect", xmin=c(0, 1.05), xmax=c(1, 2.05), ymin=rep(0.5, 2), 
           ymax=rep(0.5+nrow(d), 2), alpha=0, color="grey")+
  geom_ribbon(aes(group=group, color=as.factor(Extinct)), alpha=0, size=0.3)+
  annotate("text", label=paste0("italic(", d$Genus[prop_mu], ")"), 
           x=-0.04, y=1:nrow(d), hjust=1, vjust=0.6, parse=TRUE, size=3)+
  annotate("text", label=c("a", "b"), x=c(0.95, 2), y=rep(nrow(d), 2), size=4.5)+
  scale_x_continuous(limits=c(-0.5, 2.05))+
  scale_color_manual(values=c("blue", "red"))+
  annotate("text", label=paste0("italic('p'['i'])"),
           x=1.025, y=-0.5, vjust=1, parse=TRUE)+
  annotate("text", label=rep(seq(0, 1, by=0.25), 4), 
           x=rep(c(0, 0.25, 0.5, 0.75, 1, 1.05, 1.3, 1.55, 1.8, 2.05), 2),
           y=c(rep(nrow(d)+0.8, 10), rep(0.2, 10)), 
           vjust=(c(rep(0, 10), rep(1, 10))), 
           hjust=rep(c(0, 0.5, 0.5, 0.5, 1), 4), size=3)+
  annotate("segment", 
          x=rep(c(0, 0.25, 0.5, 0.75, 1, 1.05, 1.3, 1.55, 1.8, 2.05), 2),
          xend=rep(c(0, 0.25, 0.5, 0.75, 1, 1.05, 1.3, 1.55, 1.8, 2.05), 2),
          y=c(rep(nrow(d)+0.7, 10), rep(0.3, 10)),
          yend=c(rep(nrow(d)+0.5, 10), rep(0.5, 10)))+
  theme(legend.position="none", axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.grid=element_blank(), axis.text=element_blank(),
        axis.ticks=element_blank(), axis.title=element_blank())

ggsave("Figure_3.jpeg", plot=taxon_plot, device="jpeg", 
       dpi=1000, height=6, width=6, units="in")



#####################################################################
#####################   FIGURE 2  ###################################
#####################################################################

# Create panel for prior and posterior mu values
s <- nrow(samples_Extinct$mu)
mu_samples <- data.frame(samples=c(samples_Null$mu,
                                   samples_Extinct$mu[, 1],
                                   samples_Extinct$mu[, 2]),
                         Extinct=c(rep("Null", s),
                                   rep("Extant", s),
                                   rep("Extinct", s)))

beta_dens <- data.frame(x=seq(0, 1, length.out=200),
                        dnsty=dbeta(seq(0, 1, length.out=200),
                                    Stan_data$mu_alpha_prior,
                                    Stan_data$mu_beta_prior))
ymax <- max(sapply(unique(mu_samples$Extinct), function(x){
  dnsty <- density(mu_samples$samples[which(mu_samples$Extinct==x)], 
                   n=1e4)
  return(max(dnsty$y))
}))
panel_a <- ggplot(mu_samples, aes(x=samples))+
  annotate("line", x=beta_dens$x, y=beta_dens$dnsty, color="grey")+
  geom_density(aes(group=Extinct, color=Extinct), alpha=0.5)+
  scale_color_manual(values=c("Null"="purple", 
                              "Extant"="blue", 
                              "Extinct"="red"))+
  annotate("text", label="a", x=0.93, y=ymax*0.93, size=4.5)+
  labs(x=as.expression(bquote(italic("p"["extinct"])*","~
                              italic("p"["average"])*", and"~
                              italic("p"["extant"]))), 
       y="Density")+
  scale_x_continuous(expand=c(0.005, 0.005))+
  scale_y_continuous(expand=c(0.02, 0))+
  theme(legend.position="none", 
        panel.background=element_rect(color="black", fill="white"),
        panel.grid=element_blank())

# Create panel for prior and posterior phi values
phi_samples <- data.frame(samples=c(samples_Null$phi,
                                   samples_Extinct$phi[, 1],
                                   samples_Extinct$phi[, 2]),
                         Extinct=c(rep("Null", s),
                                   rep("Extant", s),
                                   rep("Extinct", s)))

hnorm_dens <- data.frame(x=seq(2, 2 + Stan_data$phi_scale_prior*10, 
                               length.out=1e4),
                        dnsty=dnorm(seq(2, 2 + Stan_data$phi_scale_prior*10, 
                                        length.out=1e4),
                                    2, Stan_data$phi_scale_prior),
                        mod=rep("Prior", 1e4))

phi_dens <- lapply(c("Null", "Extant", "Extinct"), function(x){
  
  df1 <- phi_samples[which(phi_samples$Extinct==x), ]
  
  df2 <- density(df1$samples, n=1e4)
  
  df <- data.frame(x=df2$x, dnsty=df2$y, mod=rep(x, length(df2$x)))
  
  return(df)
})

phi_df <- rbind(hnorm_dens, phi_dens[[1]], phi_dens[[2]], phi_dens[[3]])
phi_df <- phi_df[which(phi_df$x>=2), ]
phi_df <- phi_df[which(phi_df$dnsty > max(phi_df$dnsty)/150), ]

ymax2 <- max(phi_df$dnsty)

panel_b <- ggplot(phi_df, aes(x=x, y=dnsty, group=mod))+
  geom_line(aes(color=mod))+
  scale_color_manual(values=c("Prior"="grey", "Null"="purple", 
                              "Extant"="blue", "Extinct"="red"))+
  annotate("text", label="b", x=2+0.93*(max(phi_df$x)-min(phi_df$x)), 
           y=ymax2*0.93, size=4.5)+
  labs(x=as.expression(bquote(italic(phi["extinct"])*","~
                                italic(phi["average"])*", and"~
                                italic(phi["extant"]))), 
       y="Density")+
  scale_x_continuous(expand=c(0.005, 0.005))+
  scale_y_continuous(expand=c(0.02, 0))+
  theme(legend.position="none", 
        panel.background=element_rect(color="black", fill="white"),
        panel.grid=element_blank())

# Create panel for prior and posterior prop distribution values
sn <- 100
ns <- 300
positions <- data.frame(pos=rep(sample(1:s, sn), 3),
                        mod=c(rep("Null", sn), rep("Extant", sn), 
                              rep("Extinct", sn)))

# Get posterior mu and phi pairs
muphi_samples <- lapply(1:nrow(positions), function(x){
  if(positions$mod[x]=="Null"){
    mu <- samples_Null$mu[positions$pos[x]]
    phi <- samples_Null$phi[positions$pos[x]]
  }else if(positions$mod[x]=="Extant"){
    mu <- samples_Extinct$mu[positions$pos[x], 1]
    phi <- samples_Extinct$phi[positions$pos[x], 1]
  }else if(positions$mod[x]=="Extinct"){
    mu <- samples_Extinct$mu[positions$pos[x], 2]
    phi <- samples_Extinct$phi[positions$pos[x], 2]
  }

  return(data.frame(mod=positions$mod[x], mu=mu, phi=phi))
})
muphi_samples <- do.call("rbind", muphi_samples)


dnsty_samples <- lapply(1:nrow(muphi_samples), function(x){
  
  a <- muphi_samples$mu[x]*muphi_samples$phi[x]
  b <- (1 - muphi_samples$mu[x])*muphi_samples$phi[x]
  
  dnsty <- dbeta(seq(0, 1, length.out=ns), a, b)
  
  return(data.frame(mod=rep(muphi_samples$mod[x], ns),
                    id=as.character(rep(x, ns)),
                    x=seq(0, 1, length.out=ns), y=dnsty,
                    stringsAsFactors=FALSE))
})
dnsty_samples <- do.call("rbind", dnsty_samples)
dnsty_samples <- dnsty_samples[which(dnsty_samples$y!=Inf), ]

# Sample mu and phi values from the prior and create densities
# for plotting
dnsty_priors <- lapply(1:sn, function(x){
  
  mu <- rbeta(1, Stan_data$mu_alpha_prior, Stan_data$mu_beta_prior)
  phi <- 2 + abs(rnorm(1, 0, Stan_data$phi_scale_prior))
  
  a <- mu*phi
  b <- (1 - mu)*phi
  
  dnsty <- dbeta(seq(0, 1, length.out=ns), a, b)
  
  return(data.frame(mod=rep("Prior", ns),
                    id=as.character(rep(x+(sn*3), ns)),
                    x=seq(0, 1, length.out=ns), y=dnsty,
                    stringsAsFactors=FALSE))
})
dnsty_priors <- do.call("rbind", dnsty_priors)
dnsty_prios <- dnsty_priors[which(dnsty_priors$y!=Inf), ]
# Truncate for densities higher than those in the posteriors
dnsty_priors <- dnsty_priors[which(dnsty_priors$y < max(dnsty_samples$y)), ]


# Bind prior and posterior density samples
dnsty_df <- rbind(dnsty_priors, dnsty_samples)

# Truncate y axis for plotting
ymax <-dnsty_df[which(dnsty_df$mod!="Prior" & 
                        dnsty_df$x > 0.02 & 
                        dnsty_df$x < 0.98), ]
ymax <- max(ymax$y)
dnsty_df$y[which(dnsty_df$y > ymax)] <- ymax

# Create a panel for each posterior
for(i in 1:3){
  
  m <- c("Null", "Extant", "Extinct")[i]
  if(m=="Null"){
    lcol <- "purple"
  }else if(m=="Extant"){
    lcol <- "blue"
  }else if(m=="Extinct"){
    lcol <-"red"
  }
  
  d_df <- dnsty_df[which(dnsty_df$mod%in%c("Prior", m)), ]
  d_df$mod[which(d_df$mod==m)] <- "mod"
  
  x_title <- c("(all genera)", "(extant genera)", "(extinct genera)")[i]
  
  panel_x <- ggplot(data=NULL, aes(x=x, y=y, group=id))+
    geom_line(data=d_df[which(d_df$mod=="Prior"), ], color="grey", alpha=0.3)+
    geom_line(data=d_df[which(d_df$mod=="mod"), ], color=lcol, alpha=0.3)+
    annotate("segment", x=0, xend=1, y=ymax, yend=ymax, color="white")+
    annotate("text", label=letters[2+i], x=0.93, y=ymax*0.93, size=4.5)+
    labs(x=as.expression(bquote(italic("p"["i"])~.(x_title))), y="Density")+
    scale_x_continuous(expand=c(0.005, 0.005))+
    scale_y_continuous(expand=c(0.02, 0))+
    theme(legend.position="none", 
          panel.background=element_rect(color="black", fill="white"),
          panel.grid=element_blank())
  
  assign(paste0("panel_", letters[2+i]), panel_x)
}


# Create panel for posterior mean archaeological proportions of each
# taxon, as well as taxa sampled from the posterior.
sn <- 1000
prior_post_pos_jumble <- sample(1:(sn*2), sn*2)
yjitter_obs <- rbeta(nrow(d), 5, 5)
yjitter_post <- rbeta(sn, 5, 5)
post_pos <- sample(1:length(samples_Null$mu), sn)

prior_taxa_samples <- rbeta(sn, Stan_data$mu_alpha_prior, 
                            Stan_data$mu_beta_prior)
pts <- data.frame(dist=rep("Prior", sn),
                  x=prior_taxa_samples, y=rbeta(sn, 5, 5))

post_taxa_samples_null <- lapply(1:sn, function(y){
  x <- post_pos[y]
  a <- samples_Null$mu[x]*samples_Null$phi[x]
  b <- (1 - samples_Null$mu[x])*samples_Null$phi[x]

  df <- data.frame(dist="Null", x=rbeta(1, a, b), 
                   y=1+yjitter_post[y])
  return(df)
})
post_taxa_samples_extant <- lapply(1:sn, function(y){
  x <- post_pos[y]
  a <- samples_Extinct$mu[x, 1]*samples_Extinct$phi[x, 1]
  b <- (1 - samples_Extinct$mu[x, 1])*samples_Extinct$phi[x, 1]
  
  df <- data.frame(dist="Extant", x=rbeta(1, a, b), 
                   y=2+yjitter_post[y])
  return(df)
})
post_taxa_samples_extinct <- lapply(1:sn, function(y){
  x <- post_pos[y]
  a <- samples_Extinct$mu[x, 2]*samples_Extinct$phi[x, 2]
  b <- (1 - samples_Extinct$mu[x, 2])*samples_Extinct$phi[x, 2]
  
  df <- data.frame(dist="Extinct", x=rbeta(1, a, b), 
                   y=3+yjitter_post[y])
  return(df)
})

post_taxa_samples_null <- do.call("rbind", post_taxa_samples_null)
post_taxa_samples_extant <- do.call("rbind", post_taxa_samples_extant)
post_taxa_samples_extinct <- do.call("rbind", post_taxa_samples_extinct)
l <- list (post_taxa_samples_null, post_taxa_samples_extant,
           post_taxa_samples_extinct)

p_list <- lapply(1:3, function(x){
  p <- pts
  pts$y <- p$y + x
  df <- rbind(pts, l[[x]])
  df <- df[prior_post_pos_jumble, ]
  return(df)
})
p_df <- do.call("rbind", p_list)

# Get posterior medians for each taxon. Apply a common
# function across all posterior samples.
post_med <- function(x1, status, yaugment, mod){
  
  x <- x1$Arch_prop
  x <- x[ , which(d$Extinct%in%status)]

  df <- data.frame(Model=rep(mod, ncol(x)),
                   Extinct=as.character(d$Extinct[d$Extinct%in%status]),
                   x=apply(x, 2, median),
                   y=yaugment+yjitter_obs[which(d$Extinct%in%status)],
                   stringsAsFactors=FALSE)
}
obs_taxa <- rbind(post_med(samples_Null, c(0, 1), 1, "Null"),
                  post_med(samples_Extinct, 0, 2, "Extant"),
                  post_med(samples_Extinct, 1, 3, "Extinct"))

# Create panel f plot
panel_f <- ggplot(p_df, aes(x=x, y=y))+
  geom_point(aes(color=dist), alpha=0.5, shape=16, size=2)+
  geom_point(data=obs_taxa, aes(fill=Extinct, color=Extinct), shape=21, size=2.4)+
  scale_color_manual(values=c("Prior"="grey", "Null"="purple",
                              "Extant"="blue", "Extinct"="red",
                              "0"="white", "1"="black"))+
  scale_fill_manual(values=c("0"="black", "1"="white"))+
  annotate("text", label="f", x=0.93, y=1+(3.2*0.93), size=4.5)+
  xlab(label=as.expression(bquote(italic("p"["i"])~
                                           "(extinct, extant, and all genera)")))+
  scale_x_continuous(expand=c(0.005, 0.005))+
  scale_y_continuous(expand=c(0.02, 0), limits=c(1, 4.2))+
  theme(legend.position="none", 
        panel.background=element_rect(color="black", fill="white"),
        panel.grid=element_blank(), axis.title.y=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())

pars_plot <- plot(panel_a + panel_b + panel_c + panel_d + panel_e +
                    panel_f + plot_layout(ncol=2))

ggsave("Figure_2.jpeg", plot=pars_plot, device="jpeg", 
       dpi=1000, height=8, width=9, units="in")


#####################################################################
#####################   FIGURE 4  ###################################
#####################################################################

Null_samples <- sapply(1:length(samples_Null$mu), function(x){
  rbeta(1, samples_Null$mu[x]*samples_Null$phi[x],
        (1-samples_Null$mu[x])*samples_Null$phi[x])
})

Extinct_samples <- sapply(1:nrow(samples_Extinct$mu), function(x){
  rbeta(1, samples_Extinct$mu[x, 2]*samples_Extinct$phi[x, 2],
        (1-samples_Extinct$mu[x, 2])*samples_Extinct$phi[x, 2])
})

cat(paste0("Null model; genus >12.5 kya arch. prop. : median = ",
           round(median(Null_samples), 3), "; 95% HDI = ", 
           round(hdi(Null_samples)[[1]], 3), "-", 
           round(hdi(Null_samples)[[2]], 3), "\n"))
cat(paste0("Extinction model; genus >12.5 kya arch. prop. : median = ",
           round(median(Extinct_samples), 3), "; 95% HDI = ", 
           round(hdi(Extinct_samples)[[1]], 3), "-", 
           round(hdi(Extinct_samples)[[2]], 3), "\n"))

DF_samples <- data.frame(model=c(rep("Null", length(Null_samples)),
                                 rep("Extinct", length(Extinct_samples))),
                         x=c(Null_samples, Extinct_samples),
                         stringsAsFactors=FALSE)

Figure_4 <- ggplot(DF_samples, aes(x=x, group=model))+
  geom_density(aes(color=model), alpha=0.5)+
  scale_color_manual(values=c("Null"="purple", "Extinct"="red"))+
  labs(x=as.expression(bquote(italic("p"["i"])~
                                       "(genera last appearing before 12.5 "^14*
                                "C kyr BP)")), 
       y="Density")+
  scale_x_continuous(expand=c(0.005, 0.005), breaks=seq(0, 1, by=0.2))+
  scale_y_continuous(expand=c(0.02, 0), breaks=seq(0, 12, by=2))+
  theme(legend.position="none", 
        panel.background=element_rect(color="black", fill="white"),
        panel.grid=element_blank(),
        axis.text=element_text(size=8))

ggsave("Figure_4.jpeg", plot=Figure_4, device="jpeg", 
       dpi=1000, height=3, width=4.5, units="in")
