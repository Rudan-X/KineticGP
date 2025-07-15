rm(list = ls())
library(rBayesianOptimization)
library(DiceKriging)
library(pso)
library(Metrics)
library(reshape2)
library(ggplot2)
library(ggpubr)

path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
# path <- "/home/mpimp-golm.mpg.de/xu2004/KineticGP/"
setwd(path)

method <- "rrBLUP"
RHs <- c(0.65)# seq(0.45,0.75,0.1)
PARs <- seq(300,1500,200)

ncond <- length(RHs)*length(PARs)
topX <- 10;
scale <- "original";

all_simAsat <- matrix(nrow = 504, ncol = ncond)

conditions <- c()
for (r in 1:length(RHs)){
  filen <- paste0("results/KineticGP_Asat_env/simulatedA2223_top",topX,"_",scale,"_",method,"_RH",RHs[r],".csv")
  simAsat <- read.csv(filen)
  all_simAsat[,((r-1)*length(PARs)+1):(r*length(PARs))] <- unlist(simAsat[,3:9])
  # conditions <- c(conditions,paste0(paste0(rep("RH",length(PARs)),RHs[r]),paste0("_PAR",PARs)))
  conditions <- c(conditions,paste0("PAR",PARs))
}


ind_complete <- which(complete.cases(all_simAsat))

simulated_rates <- all_simAsat[ind_complete,]

measAQ2223<- read.csv("data/processed_data/Testing_AQcurves_years22&23_accession.csv")
measAQ2223<-measAQ2223[,1:3]

years <- measAQ2223$Year
measured_rates <- measAQ2223[ind_complete,3]
years <- years[ind_complete]

nvar <- length(measured_rates)
condition_idx=rep(1,nvar)

extract_elements <- function(index_vector, matrix) {
  # Initialize an empty vector to store the results
  result <- numeric(length(index_vector))
  
  # Loop through each element of the index_vector
  for (i in 1:length(index_vector)) {
    result[i] <- matrix[i, index_vector[i]]
  }
  
  # Return the result vector
  return(result)
}


objective_function <- function(condition_idx) {
  # Convert inputs to integer indices
  
  condition_idx <- round(condition_idx)
  
  # Calculate correlation for a given genotype and condition
  sim_values <- extract_elements(condition_idx,simulated_rates)
  meas_values <- measured_rates
  
  # Calculate correlation (you can use different measures like R-squared, RMSE, etc.)
  # obj_val <- rmse(meas_values, sim_values)
  obj_val <- - cor(sim_values, meas_values)
  
  # Return the negative of the correlation (because the optimization will minimize)
  return(obj_val)
}


result <- psoptim(par = rep(1,nvar), fn = objective_function, 
                  lower =rep(1, nvar), upper = rep(ncond, nvar), 
                  control = list(maxit = 5000, s = 100, maxit.stagnate = 500))


best_cond <- round(result$par)
field_cond <- conditions[best_cond]


# Calculate correlation for a given genotype and condition
sim_values <- extract_elements(best_cond,simulated_rates)
meas_values <- measured_rates


genotypes <- paste0("",seq(1,nvar))
dataA <- data.frame(Measured = meas_values, Simulated = sim_values, Year=years, Condition = field_cond, Genotype = genotypes, type="Adjusted PAR")
dataA$Condition <- factor(dataA$Condition, levels = paste0("PAR",PARs))

table(dataA$Condition)

g1 <- ggscatter(dataA, y="Simulated",x="Measured", color="Condition", size = 3,alpha=0.8)   +  #,label="Genotype"
  theme_bw()+
  geom_smooth(method = "lm", color = "black") +
  facet_grid(. ~ Year)+
  stat_cor(label.y=15.5,digits = 3,method = "pearson")+
  theme(legend.position = "bottom",legend.text = element_text(size=9,face="plain"),legend.title=element_text(size=10,face="bold"),
        axis.text=element_text(size=11,face="bold",color="black"), axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size=11,face="bold"))+
  labs(x="Measured", y = "Simulated")+
  scale_color_brewer(palette = "Paired")

####################################################################################

RHs <- seq(0.45,0.75,0.1) # c(0.65)
PARs <- c(300) # seq(300,1500,200)

ncond <- length(RHs)*length(PARs)
topX <- 10;
scale <- "original";

all_simAsat <- matrix(nrow = 504, ncol = ncond)

conditions <- paste0("RH",RHs)
for (r in 1:length(RHs)){
  filen <- paste0("results/KineticGP_Asat_env/simulatedA2223_top",topX,"_",scale,"_",method,"_RH",RHs[r],".csv")
  simAsat <- read.csv(filen)
  all_simAsat[,r] <- unlist(simAsat$PAR_300)
}

ind_complete <- which(complete.cases(all_simAsat))

simulated_rates <- all_simAsat[ind_complete,]
measAQ2223<- read.csv("data/processed_data/Testing_AQcurves_years22&23_accession.csv")

years <- measAQ2223$Year
measured_rates <- measAQ2223[ind_complete,3]
years <- years[ind_complete]


nvar <- length(measured_rates)
condition_idx=rep(1,nvar)

result <- psoptim(par = rep(1,nvar), fn = objective_function, 
                  lower =rep(1, nvar), upper = rep(ncond, nvar), 
                  control = list(maxit = 5000, s = 100, maxit.stagnate = 500))


best_cond <- round(result$par)
field_cond <- conditions[best_cond]

# Calculate correlation for a given genotype and condition
sim_values <- extract_elements(best_cond,simulated_rates)
meas_values <- measured_rates


genotypes <- paste0("",seq(1,nvar))
dataB <- data.frame(Measured = meas_values, Simulated = sim_values, Year=years, Condition = field_cond, Genotype = genotypes, type="Adjusted")


table(dataB$Condition)


# simAQ <- read.csv(paste0("results/KineticGP_Asimulation/simulatedA_top",topX,"_",scale,"_",method,"_field.csv"))
# simAQ <- simAQ[simAQ$Year==2021,]
# simAQ <- simAQ[ind_complete,]
# dataC <- data.frame(Measured=meas_values,Simulated=simAQ$PAR_1800,Condition="PAR300, RH0.65", Genotype = genotypes, type="Default")

g2 <- ggscatter(dataB, y="Simulated",x="Measured", color="Condition", size = 3,alpha=0.9)   +  #,label="Genotype"
  theme_bw()+
  facet_grid(. ~ Year)+
  geom_smooth(method = "lm", color = "black") +
  stat_cor(label.y=15.5,digits = 3,method = "pearson")+
  theme(legend.position = "bottom",legend.text = element_text(size=9,face="plain"),legend.title=element_text(size=10,face="bold"),
        axis.text=element_text(size=11,face="bold",color="black"), axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size=11,face="bold"))+
  labs(x="Measured", y = "Simulated")+
  scale_color_brewer(palette = "Dark2") +
  guides(color = guide_legend(nrow = 2))



####################################################################################

RHs <- seq(0.45,0.75,0.1) 
PARs <- seq(300,1500,200) 

ncond <- length(RHs)*length(PARs)
topX <- 10;
scale <- "original";

all_simAsat <- matrix(nrow = 504, ncol = ncond)

conditions <-c()
for (r in 1:length(RHs)){
  filen <- paste0("results/KineticGP_Asat_env/simulatedA2223_top",topX,"_",scale,"_",method,"_RH",RHs[r],".csv")
  simAsat <- read.csv(filen)
  all_simAsat[,((r-1)*length(PARs)+1):(r*length(PARs))] <- unlist(simAsat[,3:9])
  conditions <- c(conditions,paste0(paste0(rep("RH",length(PARs)),RHs[r]),paste0("_PAR",PARs)))
}

ind_complete <- which(complete.cases(all_simAsat))

simulated_rates <- all_simAsat[ind_complete,]


measAQ2223<- read.csv("data/processed_data/Testing_AQcurves_years22&23_accession.csv")
measAQ2223<-measAQ2223[,1:3]

years <- measAQ2223$Year
measured_rates <- measAQ2223[ind_complete,3]
years <- years[ind_complete]


nvar <- length(measured_rates)
condition_idx=rep(1,nvar)

result <- psoptim(par = rep(1,nvar), fn = objective_function, 
                  lower =rep(1, nvar), upper = rep(ncond, nvar), 
                  control = list(maxit = 5000, s = 100, maxit.stagnate = 500))


best_cond <- round(result$par)
field_cond <- conditions[best_cond]

# Calculate correlation for a given genotype and condition
sim_values <- extract_elements(best_cond,simulated_rates)
meas_values <- measured_rates


genotypes <- paste0("",seq(1,nvar))
dataC <- data.frame(Measured = meas_values, Simulated = sim_values, Year=years, Condition = field_cond, Genotype = genotypes, type="Adjusted")
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(12, "Set3"))(28)

g3 <- ggscatter(dataC, y="Simulated",x="Measured", color="Condition", size = 3)   +  #,label="Genotype"
  theme_bw() +
  facet_grid(. ~ Year) +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(label.y=15.5,digits = 3,method = "pearson")+
  theme(legend.position = "right",legend.text = element_text(size=8,face="plain"),legend.title=element_text(size=10,face="bold"),
        axis.text=element_text(size=11,face="bold",color="black"), axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size=11,face="bold"))+
  labs(x="Measured", y = "Simulated")+
  scale_color_manual(values = colors) # +
  # guides(color = guide_legend(nrow = 2))

save.image(file = "results/RData/season2223_inversion.RData")


top_row = ggarrange(g1, g2, ncol = 2, labels = c("a", "b"))
bottom_row = ggarrange( g3, labels = c("c"))
ggarrange(top_row, bottom_row , ncol = 1, heights = c(1,1) )

ggsave(paste0("Figures/SFsec6.1_inversion.png"),width =11, height = 8)
