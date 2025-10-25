############################################################
#XXXXXXXXXXXXX
#workspace prep
############################################################

#set directory
setwd(" ")

###load packages

library(performance)
library(officer)
library(flextable)
library(psych)
library(car)
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(lmtest)
library(FactoMineR)
library(lme4)  
library(nlme) 
library(factoextra)
library(vcd)
library(missMDA)
library(MuMIn)
library(survival)
library(Hmisc)
library(glmm)
library(glmmTMB)
library(Hmisc)
library(cowplot)
library(MASS)
library(lmPerm)
library(rptR)
library(GGally)
library(fitdistrplus)
library(dplyr)

#Biological data 

Data<-read.csv2("data.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)
Data_beha<-read.csv2("data_behaviors.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)


##########
#Test whether behavior shows repeatability over time in order to be considered personality.
##########
#Time_latency

plot(fitdist(Data$Time_latency,"pois"))
plot(fitdist(Data$Time_latency,"norm"))
plot(fitdist(Data$Time_latency,"nbinom"))
#Time_latency followed a negative binomial distribution
#To take into account this type of distribution, we used 
#the method of Nakagawa et al., (2017) based on the function glmer.nb 
#from the package lme4 (version 1.1-29.).


r1<-glmer.nb(Time_latency ~Session+Sex+(1|Individual)+(1|Mother), data=Data_beha)
summary(r1)
#no sex effect so we make another model without the sex factor
r2<-glmer.nb(Time_latency ~Session+(1|Individual)+(1|Mother), data=Data_beha)
summary(r2)

thetaN2 <- getME(r2, "glmer.nb.theta")
lambda2 <- as.numeric(exp(fixef(r2) + 0.5 * (as.numeric(VarCorr(r2)$Individu))))
VarOlN2 <- log(1 + (1/lambda2) + (1/thetaN2)) # log-normal approximation
VarOlN2
c(VarOlN2 = VarOlN2)
ICCrawPop2 <- as.numeric(VarCorr(r2)$Individu)/(sum(as.numeric(VarCorr(r2))) +
                                                 VarOlN2)
c(ICCrawPop2 = ICCrawPop2)

#low repeatability of 0.26 
#latency is not considered as personality traits

#########
#Number_grooming 

plot(fitdist(Data$Number_grooming ,"pois"))
plot(fitdist(Data$Number_grooming ,"norm"))
plot(fitdist(Data$Number_grooming ,"nbinom"))
#Number_grooming followed a negative binomial distribution
#To take into account this type of distribution, we used 
#the method of Nakagawa et al., (2017) based on the function glmer.nb 
#from the package lme4 (version 1.1-29.).

r3<-glmer.nb(Number_grooming ~Session+Sex+(1|Individual)+(1|Mother), data=Data_beha)
summary(r3)
#no sex ans session effect so we make another model without the sex factor 
r4<-glmer.nb(Number_grooming ~(1|Individual)+(1|Mother), data=Data_beha)
summary(r4)

thetaN4 <- getME(r4, "glmer.nb.theta")
lambda4 <- as.numeric(exp(fixef(r4) + 0.5 * (as.numeric(VarCorr(r4)$Individu))))
VarOlN4 <- log(1 + (1/lambda4) + (1/thetaN4)) # log-normal approximation
VarOlN4
c(VarOlN4 = VarOlN4)
ICCrawPop4 <- as.numeric(VarCorr(r4)$Individu)/(sum(as.numeric(VarCorr(r4))) +
                                                 VarOlN4)
c(ICCrawPop4 = ICCrawPop4)

#no repeatability
#Number_grooming is not considered as personality traits

######
#Number_transition

plot(fitdist(Data$Number_transition,"pois"))
plot(fitdist(Data$Number_transition,"norm"))
plot(fitdist(Data$Number_transition,"nbinom"))
#Number_transition followed a negative binomial distribution

r5<-glmer.nb(Number_transition ~Session+Sex+(1|Individual), data=Data_beha)
summary(r5)
#no sex and session effect so we make another model without these factors 
r6<-glmer.nb(Number_transition ~(1|Individual), data=Data_beha)
summary(r6)

thetaN6 <- getME(r6, "glmer.nb.theta")
lambda6 <- as.numeric(exp(fixef(r6) + 0.5 * (as.numeric(VarCorr(r6)$Individu))))
VarOlN6 <- log(1 + (1/lambda6) + (1/thetaN6)) # log-normal approximation
VarOlN6
c(VarOlN6 = VarOlN6)
ICCrawPop6 <- as.numeric(VarCorr(r6)$Individu)/(sum(as.numeric(VarCorr(r6))) +
                                                 VarOlN6)
c(ICCrawPop6 = ICCrawPop6)

#moderate repeatability of 0.37 
#Number_transition is considered as personality traits

########
#Number_rearing

plot(fitdist(Data$Number_rearing,"pois"))
plot(fitdist(Data$Number_rearing,"norm"))
plot(fitdist(Data$Number_rearing,"nbinom"))
#Number_rearing followed a negative binomial distribution

r7<-glmer.nb(Number_rearing ~Session+Sex+(1|Individual)+(1|Mother), data=Data_beha)
summary(r7)
#no sex and session effect so we make another model without these factors 
r8<-glmer.nb(Number_rearing ~(1|Individual)+(1|Mother), data=Data_beha)
summary(r8)

thetaN8 <- getME(r8, "glmer.nb.theta")
lambda8 <- as.numeric(exp(fixef(r8) + 0.5 * (as.numeric(VarCorr(r8)$Individu))))
VarOlN8 <- log(1 + (1/lambda8) + (1/thetaN8)) # log-normal approximation
VarOlN8
c(VarOlN8 = VarOlN8)
ICCrawPop8 <- as.numeric(VarCorr(r8)$Individu)/(sum(as.numeric(VarCorr(r8))) +
                                                 VarOlN8)
c(ICCrawPop8 = ICCrawPop8)

#no repeatability
#Number_rearing is not considered as personality traits


# ============================================================
# === 1. Define variables and categories ====================
# ============================================================

traits <- c(
  "Growth_rate", "Temperature_torpor", "Torpor_bout_duration",
  "Time_inter_torpor", "Time_torpor", "Delta_telomere_experiment",
  "Telomere_pre_hibernation",
  "log_Mean_cortisol", "Number_transition",
  "Offspring_per_litter", "Offspring_number", "Offspring_growth_rate"
)


# Create a vector to store the names of the transformed variables
variables <- c()

for (var in traits) {
  # Subset data without missing values
  temp_data <- na.omit(Data[, c(var, "Mother", "Sex")])
  
  # Fit the LME model with Sex as a fixed effect
  model <- lme(as.formula(paste(var, "~ Sex")), random = ~1|Mother, method = "ML", data = temp_data)
  
  # Store the residuals in the main Data frame
  res_name <- paste0("res_", var)
  Data[, res_name] <- NA
  Data[rownames(temp_data), res_name] <- residuals(model)
  
  # Add the name of the transformed variable to 'variables'
  variables <- c(variables, res_name)
  
  cat(sprintf("Residuals for %s computed: %d values\n", var, sum(!is.na(residuals(model)))))
}



categories <- list(
  Reproduction = c("res_Offspring_per_litter", "res_Offspring_growth_rate", "res_Offspring_number"),
  Hibernation = c( "res_Temperature_torpor","res_Torpor_bout_duration", 
                   "res_Time_inter_torpor", "res_Time_torpor"),
  Telomere = c("res_Delta_telomere_experiment", "res_Telomere_pre_hibernation"),
  Cortisol = c("res_log_Mean_cortisol"),
  Transition = c("res_Number_transition"),
  Growth = c("res_Growth_rate")
)
# ============================================================
# === 2. Utility functions ==================================
# ============================================================
calculate_kmo <- function(data, vars) {
  kmo_result <- KMO(data[, vars, drop = FALSE])
  return(kmo_result$MSA[1])
}

test_bartlett <- function(data, vars) {
  corr_mat <- cor(data[, vars, drop = FALSE])
  if(det(corr_mat) < 1e-10) return(1)
  bartlett_test <- cortest.bartlett(corr_mat, n = nrow(data))
  return(bartlett_test$p.value)
}

check_categories <- function(vars, categories) {
  used_cats <- sapply(vars, function(v) {
    cat_name <- names(Filter(function(x) v %in% x, categories))
    if(length(cat_name)==0) return(NA)
    return(cat_name)
  })
  return(length(unique(used_cats)) == length(vars))
}

# ============================================================
# === 3. Generate all valid combinations ====================
# ============================================================
vars_to_check <- intersect(variables, names(Data))
combinations <- unlist(
  lapply(2:length(vars_to_check), function(i) combn(vars_to_check, i, simplify=FALSE)),
  recursive=FALSE
)

valid_combinations <- Filter(function(x) check_categories(x, categories), combinations)

kmo_scores <- c()
valid_sets  <- list()

# ============================================================
# === 4. Optimize variable selection (KMO + Bartlett) =======
# ============================================================
for (comb in valid_combinations) {
  subset <- na.omit(Data[, comb, drop=FALSE])
  
  bartlett_p <- tryCatch(test_bartlett(subset, comb), error=function(e) NA)
  if(is.na(bartlett_p) || bartlett_p >= 0.05) next
  
  kmo_val <- tryCatch(calculate_kmo(subset, comb), error=function(e) NA)
  if(!is.na(kmo_val)) {
    kmo_scores <- c(kmo_scores, kmo_val)
    valid_sets <- append(valid_sets, list(comb))
  }
}

# ============================================================
# === 5. Identify best combination ==========================
# ============================================================
if(length(kmo_scores) > 0) {
  best_index <- which.max(kmo_scores)
  best_combination <- valid_sets[[best_index]]
  best_kmo <- max(kmo_scores)
  
  cat("✅ Best variable combination:", paste(best_combination, collapse=", "), "\n")
  cat("✅ Best KMO score:", round(best_kmo,4), "\n")
} else {
  stop("No valid combination according to Bartlett’s test.")
}

# ============================================================
# === 6 PC1 extraction =================
# ============================================================
vars_pca <- best_combination
data_pca <- na.omit(Data[, vars_pca, drop=FALSE])
data_scaled <- scale(data_pca)

res_pca <- PCA(data_scaled, scale.unit=FALSE, ncp=2, graph=FALSE)
Data$Syndrome_PC1 <- NA
Data[rownames(data_pca), "Syndrome_PC1"] <- res_pca$ind$coord[,1]
res_pca_syndrome$var$coord
rownames(res_pca_syndrome$var$coord) <- c("Growth rate", "Offspring growth rate", "Number of transitions") 

fviz_pca_var(
   res_pca_syndrome, 
  axes = c(1, 2),
  col.var = "red",
  repel = TRUE,
  arrowsize = 1.2,        
  labelsize = 6,          
  title = "Biplot PCA: variables only",
  circle = TRUE
) + 
  coord_equal() +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_line(color = "grey90"),
    axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.margin = margin(10, 10, 10, 10)
  )


# ============================================================
# === 7. PCA robustness check via permutation test ==========
# ============================================================
n_perm <- 500
n_pc <- 2
perm_var <- matrix(NA, nrow=n_perm, ncol=n_pc)
var_obs <- res_pca$eig[1:n_pc,2]

for(p in 1:n_perm) {
  perm_data <- apply(data_scaled, 2, sample)
  res_perm <- tryCatch(PCA(perm_data, scale.unit=FALSE, ncp=n_pc, graph=FALSE), error=function(e) NULL)
  if(!is.null(res_perm)) perm_var[p,] <- res_perm$eig[1:n_pc,2]
}

p_values <- sapply(1:n_pc, function(k) mean(perm_var[,k] >= var_obs[k], na.rm=TRUE))

for(k in 1:n_pc) {
  cat(sprintf("PC%d: Observed variance = %.2f%%, Permutation p-value = %.3f\n", k, var_obs[k], p_values[k]))
}

# ============================================================
# === 8. LME for PC1 ~ Littersize, Birth, Sex =============
# ============================================================
Datana <- na.omit(Data[, c("Syndrome_PC1","Littersize","Birth","Mother","Sex")])

model_full <- lme(Syndrome_PC1 ~ Littersize * Sex + Birth * Sex,
                  random = ~1 | Mother, data=Datana, method="ML")

# Model selection using dredge (AICc)
model_selection <- dredge(model_full, fixed = ~ +(1|Mother), rank="AICc", m.lim=c(NA,4))
summary(model.avg(model_selection, delta<5))

# Simplified model
model_simple <- lme(Syndrome_PC1 ~ Littersize * Sex,
                    random = ~1 | Mother, data=Datana, method="ML")
summary(model_simple)

# ============================================================
# === 9. Model diagnostics ===================================
# ============================================================
residuals <- resid(model_simple)
plot(density(residuals), main="Residual Density", xlab="Residuals")
curve(dnorm(x, mean=mean(residuals), sd=sd(residuals)), col="red", lwd=2, add=TRUE)
bptest(model_simple)
performance::r2(model_simple)

# ============================================================
# === 10. Visualization =====================================
# ============================================================
ggplot(Datana, aes(x=Littersize, y=Syndrome_PC1, color=Sex)) +
  geom_point(size=3, alpha=0.8) +
  geom_smooth(method="lm", se=TRUE, aes(fill=Sex)) +
  scale_color_manual(values=c("F"="#E41A1C", "M"="#377EB8")) +
  scale_fill_manual(values=c("F"="#E41A1C", "M"="#377EB8")) +
  theme_classic(base_size=14) +
  labs(x="Litter size", y="Syndrome PC1", color="Sex", fill="Sex") +
  theme(legend.position="top")



# ============================================================
# === Female-only PCA and LME =================================
# ============================================================
DataF <- subset(Data, Sex == "F")
variables_F <- c("Growth_rate", "Delta_telomere_experiment", "Telomere_pre_hibernation",
               "Temperature_torpor", "Torpor_bout_duration", "Time_inter_torpor", "Time_torpor",
               "log_Mean_cortisol", "Offspring_per_litter", "Offspring_growth_rate",
               "Offspring_number", "Number_transition")

categories_F <- list(
  Reproduction = c("Offspring_per_litter", "Offspring_growth_rate", "Offspring_number"),
  Hibernation = c("Temperature_torpor","Torpor_bout_duration","Time_inter_torpor","Time_torpor"),
  Telomere = c("Delta_telomere_experiment","Telomere_pre_hibernation"),
  Cortisol = c("log_Mean_cortisol"),
  Transition = c("Number_transition"),
  Growth = c("Growth_rate")
)
# Generate all valid combinations
vars_to_check <- intersect(variables_F, names(DataF))
combinations <- unlist(
  lapply(2:length(vars_to_check), function(i) combn(vars_to_check, i, simplify = FALSE)),
  recursive = FALSE
)
valid_combinations <- Filter(function(x) check_categories(x, categories_F), combinations)

kmo_scores <- c()
best_sets <- list()

for (comb in valid_combinations) {
  subset <- na.omit(DataF[, comb, drop=FALSE])
  bartlett_p <- tryCatch(test_bartlett(subset, comb), error=function(e) NA)
  if (is.na(bartlett_p) || bartlett_p >= 0.05) next
  kmo_val <- tryCatch(calculate_kmo(subset, comb), error=function(e) NA)
  if (!is.na(kmo_val)) {
    kmo_scores <- c(kmo_scores, kmo_val)
    best_sets <- append(best_sets, list(comb))
  }
}

# Select best combination
best_index <- which.max(kmo_scores)
best_combination <- best_sets[[best_index]]
best_kmo <- max(kmo_scores)

cat("✅ Female best variables:", paste(best_combination, collapse=", "), "\n")
cat("✅ KMO:", round(best_kmo,4), "\n")

# PCA 
data_pc1 <- na.omit(DataF[, best_combination, drop=FALSE])
data_scaled <- scale(data_pc1)
res_pcaF <- PCA(data_scaled, scale.unit=FALSE, ncp=2, graph=FALSE)
DataF$Syndrome_PC1 <- NA
DataF[rownames(data_pc1), "Syndrome_PC1"] <- res_pcaF$ind$coord[,1]

res_pcaF$var$coord
rownames(res_pcaF$var$coord) <- c("Growth rate", "Offspring growth rate", "Number of transitions") 
fviz_pca_var(
  res_pcaF, 
  axes = c(1, 2),
  col.var = "red",
  repel = TRUE,
  arrowsize = 1.2,        
  labelsize = 6,          
  title = "Biplot PCA: variables only",
  circle = TRUE
) + 
  coord_equal() +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_line(color = "grey90"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.margin = margin(10, 10, 10, 10)
  )



# Permutation test
n_perm <- 500
n_pc <- 2
perm_var <- matrix(NA, nrow=n_perm, ncol=n_pc)
var_obs <- res_pcaF$eig[1:n_pc,2]

for (p in 1:n_perm) {
  perm_data <- apply(data_scaled, 2, sample)
  res_perm <- tryCatch(PCA(perm_data, scale.unit=FALSE, ncp=n_pc, graph=FALSE), error=function(e) NULL)
  if(!is.null(res_perm)) perm_var[p,] <- res_perm$eig[1:n_pc,2]
}

p_values <- sapply(1:n_pc, function(k) mean(perm_var[,k] >= var_obs[k], na.rm=TRUE))
for(k in 1:n_pc) cat(sprintf("Female PC%d: Observed variance = %.2f%%, Perm p = %.3f\n", k, var_obs[k], p_values[k]))

# LME for females
DataFna_model <- na.omit(DataF[, c("Syndrome_PC1", "Littersize", "Birth", "Mother")])

# Full model
model_full <- lme(Syndrome_PC1 ~ Littersize + Birth, random = ~1|Mother, method="ML", data=DataFna_model)

# Optional: model selection
summary(model.avg(dredge(model_full, fixed=~+(1|Mother), rank="AICc", m.lim=c(NA,4)), delta<5))

# Simple model
model_simple <- lme(Syndrome_PC1 ~ Littersize, random=~1|Mother, method="ML", data=DataFna_model)
summary(model_simple)

#Model diagnostics ---

# Extract residuals
residuals <- resid(model_simple)

# Check residual normality
plot(density(residuals), main = "Residual Density", xlab = "Residuals")
curve(dnorm(x, mean = mean(residuals), sd = sd(residuals)), col = "red", lwd = 2, add = TRUE)

# Test homoscedasticity
bptest(model_simple)

# Compute R²
performance::r2(model_simple)


# --------------------------------------
# 9. Visualization
# --------------------------------------

ggplot(DataFna_model, aes(x = Syndrome_PC1, y = Littersize)) +
  geom_jitter(width = 0.05, height = 0.05, size = 4, alpha = 0.6, color = "steelblue") +  # jitter
  geom_smooth(method = "lm", se = TRUE, color = "darkred", size = 1.2) +
  theme_minimal(base_size = 14) +
  labs(
    x = "PC1F",
    y = "Litter size",
    title = "Relationship between PC1 and Litter Size"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )


DataFna_model
# ============================================================
# === Male-only PCA and LME ==================================
# ============================================================
DataM <- subset(Data, Sex == "M")


# Variables of interest
variables_M <- c(
  "Growth_rate", "Delta_telomere_experiment", "Telomere_pre_hibernation",
  "Temperature_torpor", "log_Mean_cortisol", "Torpor_bout_duration",
  "Time_inter_torpor", "Time_torpor", "Offspring_per_litter",
  "Offspring_growth_rate", "Offspring_number", "Number_transition",
  "Testosterone_post_hibernation"
)

# Conceptual categories
categories_M <- list(
  Reproduction = c("Offspring_per_litter", "Offspring_growth_rate", "Offspring_number", "Testosterone_post_hibernation"),
  Hibernation  = c("Time_torpor", "Torpor_bout_duration", "Time_inter_torpor", "Temperature_torpor"),
  Telomere     = c("Delta_telomere_experiment", "Telomere_pre_hibernation"),
  Cortisol     = c("log_Mean_cortisol"),
  Transition   = c("Number_transition"),
  Growth       = c("Growth_rate")
)


# Generate all valid combinations (2–3 variables, at least 1 Reproduction)
has_reproduction <- function(vars, categories) any(vars %in% categories$Reproduction)
vars_to_check <- intersect(variables_M, names(DataM))
combinations <- unlist(lapply(2:3, function(i) combn(vars_to_check, i, simplify=FALSE)), recursive=FALSE)
valid_combinations <- Filter(function(x) check_categories(x, categories_M) && has_reproduction(x, categories_M), combinations)

kmo_scores <- c()
best_sets <- list()
for(comb in valid_combinations) {
  subset <- na.omit(DataM[, comb, drop=FALSE])
  if(nrow(subset) < 10 || ncol(subset) < 2) next
  bartlett_p <- tryCatch(test_bartlett(subset, comb), error=function(e) NA)
  if(is.na(bartlett_p) || bartlett_p >= 0.05) next
  kmo_val <- tryCatch(calculate_kmo(subset, comb), error=function(e) NA)
  if(!is.na(kmo_val)) {
    kmo_scores <- c(kmo_scores, kmo_val)
    best_sets <- append(best_sets, list(comb))
  }
}

best_index <- which.max(kmo_scores)
best_combination <- best_sets[[best_index]]
best_kmo <- max(kmo_scores)
cat("✅ Male best variables:", paste(best_combination, collapse=", "), "\n")
cat("✅ KMO:", round(best_kmo,4), "\n")

# PCA 
data_pc1 <- na.omit(DataM[, best_combination, drop=FALSE])
data_scaled <- scale(data_pc1)
res_pcaM <- PCA(data_scaled, scale.unit=FALSE, ncp=2, graph=FALSE)
DataM$Syndrome_PC1 <- NA
DataM[rownames(data_pc1), "Syndrome_PC1"] <- res_pcaM$ind$coord[,1]
res_pcaM$var$coord
rownames(res_pcaM$var$coord) <- c("Telomere length variation", "log(Mean cortisol)", "Offspring growth rate") 
fviz_pca_var(
  res_pcaM, 
  axes = c(1, 2),
  col.var = "red",
  repel = TRUE,
  arrowsize = 1.2,        
  labelsize = 6,          
  title = "Biplot PCA: variables only",
  circle = TRUE
) + 
  coord_equal() +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_line(color = "grey90"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.margin = margin(10, 10, 10, 10)
  )


# Permutation test
perm_var <- matrix(NA, nrow=n_perm, ncol=n_pc)
var_obs <- res_pcaM$eig[1:n_pc,2]
for(p in 1:n_perm) {
  perm_data <- apply(data_scaled, 2, sample)
  res_perm <- tryCatch(PCA(perm_data, scale.unit=FALSE, ncp=n_pc, graph=FALSE), error=function(e) NULL)
  if(!is.null(res_perm)) perm_var[p,] <- res_perm$eig[1:n_pc,2]
}
p_values <- sapply(1:n_pc, function(k) mean(perm_var[,k] >= var_obs[k], na.rm=TRUE))
for(k in 1:n_pc) cat(sprintf("Male PC%d: Observed variance = %.2f%%, Perm p = %.3f\n", k, var_obs[k], p_values[k]))

#No PCA valide for males



# ============================================================
# ============================================================
# ============================================================
# === Same test with hibernation pattern correction ==========
# ============================================================
# ============================================================
# ============================================================



# --- 1. Define traits ---
# ============================================================
# Compute residuals and store variable names in 'variables'
# ============================================================

# ============================================================
# Compute residuals for all traits and store transformed names
# ============================================================
Data<-read.csv2("data.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)


all_traits <- c(
  "Growth_rate", "Delta_telomere_experiment", "Telomere_pre_hibernation",
  "log_Mean_cortisol", "Number_transition",
  "Offspring_per_litter", "Offspring_number", "Offspring_growth_rate",
  "Temperature_torpor", "Torpor_bout_duration", "Time_inter_torpor", "Time_torpor"
)

hibernation_traits <- c("Temperature_torpor","Torpor_bout_duration","Time_inter_torpor","Time_torpor")

# Create empty vector to store names of transformed variables
variables <- c()

for (var in all_traits) {
  
  # Determine fixed effects formula and columns to keep
  if (var %in% hibernation_traits) {
    fixed_effects <- "~ Sex + Hamster_location"
    temp_data <- na.omit(Data[, c(var, "Mother", "Sex", "Hamster_location")])
  } else {
    fixed_effects <- "~ Sex"
    temp_data <- na.omit(Data[, c(var, "Mother", "Sex")])
  }
  
  # Fit LME
  model <- lme(as.formula(paste(var, fixed_effects)),
               random = ~1 | Mother,
               method = "ML",
               data = temp_data)
  
  # Store residuals in Data
  res_name <- paste0("res_", var)
  Data[, res_name] <- NA
  Data[rownames(temp_data), res_name] <- residuals(model)
  
  # Append transformed variable name to variables vector
  variables <- c(variables, res_name)
  
  cat(sprintf("Residuals for %s computed: %d values\n", var, sum(!is.na(residuals(model)))))
}

# Check resulting variables


categories <- list(
  Reproduction = c("res_Offspring_per_litter", "res_Offspring_growth_rate", "res_Offspring_number"),
  Hibernation = c( "res_Temperature_torpor","res_Torpor_bout_duration", 
                   "res_Time_inter_torpor", "res_Time_torpor"),
  Telomere = c("res_Delta_telomere_experiment", "res_Telomere_pre_hibernation"),
  Cortisol = c("res_log_Mean_cortisol"),
  Transition = c("res_Number_transition"),
  Growth = c("res_Growth_rate")
)

# ============================================================
# === 2. Utility functions ==================================
# ============================================================
calculate_kmo <- function(data, vars) {
  kmo_result <- KMO(data[, vars, drop = FALSE])
  return(kmo_result$MSA[1])
}

test_bartlett <- function(data, vars) {
  corr_mat <- cor(data[, vars, drop = FALSE])
  if(det(corr_mat) < 1e-10) return(1)
  bartlett_test <- cortest.bartlett(corr_mat, n = nrow(data))
  return(bartlett_test$p.value)
}

check_categories <- function(vars, categories) {
  used_cats <- sapply(vars, function(v) {
    cat_name <- names(Filter(function(x) v %in% x, categories))
    if(length(cat_name)==0) return(NA)
    return(cat_name)
  })
  return(length(unique(used_cats)) == length(vars))
}

# ============================================================
# === 3. Generate all valid combinations ====================
# ============================================================
vars_to_check <- intersect(variables, names(Data))
combinations <- unlist(
  lapply(2:length(vars_to_check), function(i) combn(vars_to_check, i, simplify=FALSE)),
  recursive=FALSE
)

valid_combinations <- Filter(function(x) check_categories(x, categories), combinations)

kmo_scores <- c()
valid_sets  <- list()

# ============================================================
# === 4. Optimize variable selection (KMO + Bartlett) =======
# ============================================================
for (comb in valid_combinations) {
  subset <- na.omit(Data[, comb, drop=FALSE])
  
  bartlett_p <- tryCatch(test_bartlett(subset, comb), error=function(e) NA)
  if(is.na(bartlett_p) || bartlett_p >= 0.05) next
  
  kmo_val <- tryCatch(calculate_kmo(subset, comb), error=function(e) NA)
  if(!is.na(kmo_val)) {
    kmo_scores <- c(kmo_scores, kmo_val)
    valid_sets <- append(valid_sets, list(comb))
  }
}

# ============================================================
# === 5. Identify best combination ==========================
# ============================================================
if(length(kmo_scores) > 0) {
  best_index <- which.max(kmo_scores)
  best_combination <- valid_sets[[best_index]]
  best_kmo <- max(kmo_scores)
  
  cat("✅ Best variable combination:", paste(best_combination, collapse=", "), "\n")
  cat("✅ Best KMO score:", round(best_kmo,4), "\n")
} else {
  stop("No valid combination according to Bartlett’s test.")
}
#The PCA is the same between corrected and uncorrected hibernation variables.



# ============================================================
# === Female-only PCA and LME =================================
# ============================================================

# ============================================================
# 1. Corrected hibernation data 
# ============================================================
Data<-read.csv2("data.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)

DataF <- subset(Data, Sex == "F")


hibernation_traits <- c("Temperature_torpor","Torpor_bout_duration","Time_inter_torpor","Time_torpor")
res_hibernation <- c()  
for (var in hibernation_traits) {
 
  temp_data <- na.omit(DataF[, c(var, "Mother", "Hamster_location")])
  
  model <- lme(as.formula(paste(var, "~ Hamster_location")),
               random = ~1 | Mother,
               method = "ML",
               data = temp_data)
  
  res_name <- paste0("res_", var)
  DataF[, res_name] <- NA
  DataF[rownames(temp_data), res_name] <- residuals(model)
  
  res_hibernation <- c(res_hibernation, res_name)
  
  cat(sprintf("Residuals for %s computed: %d values\n", var, sum(!is.na(residuals(model)))))
}

variables_F <- c(
  "Growth_rate", "Delta_telomere_experiment", "Telomere_pre_hibernation",
  "res_Temperature_torpor", "res_Torpor_bout_duration", "res_Time_inter_torpor", "res_Time_torpor",
  "log_Mean_cortisol", "Offspring_per_litter", "Offspring_growth_rate",
  "Offspring_number", "Number_transition"
)


variables_F <- c(
  "Growth_rate",
  "Delta_telomere_experiment",
  "Telomere_pre_hibernation",
  "res_Temperature_torpor", 
  "res_Torpor_bout_duration", 
  "res_Time_inter_torpor", 
  "res_Time_torpor",
  "log_Mean_cortisol",
  "Offspring_per_litter",
  "Offspring_growth_rate",
  "Offspring_number",
  "Number_transition"
)

# Generate all valid combinations
vars_to_check <- intersect(variables_F, names(DataF))
combinations <- unlist(
  lapply(2:length(vars_to_check), function(i) combn(vars_to_check, i, simplify = FALSE)),
  recursive = FALSE
)
valid_combinations <- Filter(function(x) check_categories(x, categories_F), combinations)

kmo_scores <- c()
best_sets <- list()

for (comb in valid_combinations) {
  subset <- na.omit(DataF[, comb, drop=FALSE])
  bartlett_p <- tryCatch(test_bartlett(subset, comb), error=function(e) NA)
  if (is.na(bartlett_p) || bartlett_p >= 0.05) next
  kmo_val <- tryCatch(calculate_kmo(subset, comb), error=function(e) NA)
  if (!is.na(kmo_val)) {
    kmo_scores <- c(kmo_scores, kmo_val)
    best_sets <- append(best_sets, list(comb))
  }
}

# Select best combination
best_index <- which.max(kmo_scores)
best_combination <- best_sets[[best_index]]
best_kmo <- max(kmo_scores)

cat("✅ Female best variables:", paste(best_combination, collapse=", "), "\n")
cat("✅ KMO:", round(best_kmo,4), "\n")

#The PCA is the same between corrected and uncorrected hibernation variables.

# ============================================================
# === Male-only PCA and LME ==================================
# ============================================================
Data<-read.csv2("data.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)

DataM <- subset(Data, Sex == "M")


hibernation_traits <- c("Temperature_torpor","Torpor_bout_duration","Time_inter_torpor","Time_torpor")
res_hibernation <- c()  

for (var in hibernation_traits) {
  temp_data <- na.omit(DataM[, c(var, "Mother", "Hamster_location")])
  
  model <- lme(as.formula(paste(var, "~ Hamster_location")),
               random = ~1 | Mother,
               method = "ML",
               data = temp_data)
  
  res_name <- paste0("res_", var)
  DataM[, res_name] <- NA
  DataM[rownames(temp_data), res_name] <- residuals(model)
  
  res_hibernation <- c(res_hibernation, res_name)
  
  cat(sprintf("Residuals for %s computed: %d values\n", var, sum(!is.na(residuals(model)))))
}

variables_M <- c(
  "Growth_rate", "Delta_telomere_experiment", "Telomere_pre_hibernation",
  "res_Temperature_torpor", "res_Torpor_bout_duration", "res_Time_inter_torpor", "res_Time_torpor",
  "log_Mean_cortisol", "Offspring_per_litter", "Offspring_growth_rate",
  "Offspring_number", "Number_transition"
)


variables_M <- c(
  "Growth_rate",
  "Delta_telomere_experiment",
  "Telomere_pre_hibernation",
  "res_Temperature_torpor", 
  "res_Torpor_bout_duration", 
  "res_Time_inter_torpor", 
  "res_Time_torpor",
  "log_Mean_cortisol",
  "Offspring_per_litter",
  "Offspring_growth_rate",
  "Offspring_number",
  "Number_transition"
)

# Generate all valid combinations (2–3 variables, at least 1 Reproduction)
has_reproduction <- function(vars, categories) any(vars %in% categories$Reproduction)
vars_to_check <- intersect(variables_M, names(DataM))
combinations <- unlist(lapply(2:3, function(i) combn(vars_to_check, i, simplify=FALSE)), recursive=FALSE)
valid_combinations <- Filter(function(x) check_categories(x, categories_M) && has_reproduction(x, categories_M), combinations)

kmo_scores <- c()
best_sets <- list()
for(comb in valid_combinations) {
  subset <- na.omit(DataM[, comb, drop=FALSE])
  if(nrow(subset) < 10 || ncol(subset) < 2) next
  bartlett_p <- tryCatch(test_bartlett(subset, comb), error=function(e) NA)
  if(is.na(bartlett_p) || bartlett_p >= 0.05) next
  kmo_val <- tryCatch(calculate_kmo(subset, comb), error=function(e) NA)
  if(!is.na(kmo_val)) {
    kmo_scores <- c(kmo_scores, kmo_val)
    best_sets <- append(best_sets, list(comb))
  }
}

best_index <- which.max(kmo_scores)
best_combination <- best_sets[[best_index]]
best_kmo <- max(kmo_scores)
cat("✅ Male best variables:", paste(best_combination, collapse=", "), "\n")
cat("✅ KMO:", round(best_kmo,4), "\n")


#The PCA is the same between corrected and uncorrected hibernation variables.


# ============================================================
# ============================================================
# ============================================================
# === 1. Test without the constraint of ACP with reproduction 
# ============================================================
# ============================================================
# ============================================================

# ============================================================
# === 1. Load data and compute residuals ====================
# ============================================================

Data<-read.csv2("data.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)

traits <- c(
  "Growth_rate", "Temperature_torpor", "Torpor_bout_duration",
  "Time_inter_torpor", "Time_torpor", "Delta_telomere_experiment",
  "Telomere_pre_hibernation",
  "log_Mean_cortisol", "Number_transition",
  "Offspring_per_litter", "Offspring_number", "Offspring_growth_rate"
)

variables <- c()

for (var in traits) {
  temp_data <- na.omit(Data[, c(var, "Mother", "Sex")])
  model <- lme(as.formula(paste(var, "~ Sex")), random = ~1|Mother, method="ML", data=temp_data)
  
  res_name <- paste0("res_", var)
  Data[, res_name] <- NA
  Data[rownames(temp_data), res_name] <- residuals(model)
  
  variables <- c(variables, res_name)
  cat(sprintf("Residuals for %s computed: %d values\n", var, sum(!is.na(residuals(model)))))
}

categories <- list(
  Reproduction = c("res_Offspring_per_litter", "res_Offspring_growth_rate", "res_Offspring_number"),
  Hibernation  = c("res_Temperature_torpor", "res_Torpor_bout_duration", "res_Time_inter_torpor", "res_Time_torpor"),
  Telomere     = c("res_Delta_telomere_experiment", "res_Telomere_pre_hibernation"),
  Cortisol     = c("res_log_Mean_cortisol"),
  Transition   = c("res_Number_transition"),
  Growth       = c("res_Growth_rate")
)

# ============================================================
# === 2. Utility functions ==================================
# ============================================================
calculate_kmo <- function(data, vars) {
  kmo_result <- KMO(data[, vars, drop=FALSE])
  return(kmo_result$MSA[1])
}

test_bartlett <- function(data, vars) {
  corr_mat <- cor(data[, vars, drop=FALSE])
  if(det(corr_mat) < 1e-10) return(1)
  bartlett_test <- cortest.bartlett(corr_mat, n=nrow(data))
  return(bartlett_test$p.value)
}

check_categories <- function(vars, categories) {
  used_cats <- sapply(vars, function(v) {
    cat_name <- names(Filter(function(x) v %in% x, categories))
    if(length(cat_name) == 0) return(NA)
    return(cat_name)
  })
  return(length(unique(used_cats)) == length(vars))
}

# ============================================================
# === 3. Generate all valid combinations ====================
# ============================================================
vars_to_check <- intersect(variables, names(Data))
combinations <- unlist(
  lapply(2:length(vars_to_check), function(i) combn(vars_to_check, i, simplify=FALSE)),
  recursive=FALSE
)

# Keep only combinations respecting one variable per category
valid_combinations <- Filter(function(x) check_categories(x, categories), combinations)

# ============================================================
# === 4. Compute KMO and Bartlett ===========================
# ============================================================
kmo_scores <- c()
valid_sets <- list()

for(comb in valid_combinations){
  subset <- na.omit(Data[, comb, drop=FALSE])
  
  # Bartlett test first
  bartlett_p <- tryCatch(test_bartlett(subset, comb), error=function(e) NA)
  if(is.na(bartlett_p) || bartlett_p >= 0.05) next  # skip unsuitable correlation matrices
  
  # KMO calculation
  kmo_val <- tryCatch(calculate_kmo(subset, comb), error=function(e) NA)
  if(!is.na(kmo_val)){
    kmo_scores <- c(kmo_scores, kmo_val)
    valid_sets <- append(valid_sets, list(comb))
  }
}

# ============================================================
# === 5. Identify best combination ==========================
# ============================================================
if(length(kmo_scores) > 0){
  best_index <- which.max(kmo_scores)
  best_combination <- valid_sets[[best_index]]
  best_kmo <- max(kmo_scores)
  
  cat("✅ Best variable combination:", paste(best_combination, collapse=", "), "\n")
  cat("✅ Best KMO score:", round(best_kmo, 4), "\n")
} else {
  stop("No valid combination according to Bartlett’s test and category constraint.")
}

#The PCA is the same with or without the constraint of reproduction


# ============================================================
# === Female-only PCA and LME =================================
# ============================================================
# ============================================================
# === 1. Load data and define variables =====================
# ============================================================

Data<-read.csv2("data.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)

DataF <- subset(Data, Sex == "F")

variables_F <- c(
  "Growth_rate", "Delta_telomere_experiment", "Telomere_pre_hibernation",
  "Temperature_torpor", "Torpor_bout_duration", "Time_inter_torpor", "Time_torpor",
  "log_Mean_cortisol", "Offspring_per_litter", "Offspring_growth_rate",
  "Offspring_number", "Number_transition"
)

categories_F <- list(
  Reproduction = c("Offspring_per_litter", "Offspring_growth_rate", "Offspring_number"),
  Hibernation  = c("Temperature_torpor", "Torpor_bout_duration", "Time_inter_torpor", "Time_torpor"),
  Telomere     = c("Delta_telomere_experiment", "Telomere_pre_hibernation"),
  Cortisol     = c("log_Mean_cortisol"),
  Transition   = c("Number_transition"),
  Growth       = c("Growth_rate")
)

# ============================================================
# === 2. Utility functions ==================================
# ============================================================
calculate_kmo <- function(data, vars) {
  kmo_result <- KMO(data[, vars, drop=FALSE])
  return(kmo_result$MSA[1])
}

test_bartlett <- function(data, vars) {
  corr_mat <- cor(data[, vars, drop=FALSE])
  if (det(corr_mat) < 1e-10) return(1)  # skip if determinant too small (singular matrix)
  bartlett_test <- cortest.bartlett(corr_mat, n=nrow(data))
  return(bartlett_test$p.value)
}

check_categories <- function(vars, categories) {
  used_cats <- sapply(vars, function(v) {
    cat_name <- names(Filter(function(x) v %in% x, categories))
    if (length(cat_name) == 0) return(NA)
    return(cat_name)
  })
  return(length(unique(used_cats)) == length(vars))  # TRUE if all vars from different categories
}

# ============================================================
# === 3. Generate all valid combinations ====================
# ============================================================
vars_to_check <- intersect(variables_F, names(DataF))

# Generate all combinations of at least 2 variables
combinations <- unlist(
  lapply(2:length(vars_to_check), function(i) combn(vars_to_check, i, simplify=FALSE)),
  recursive=FALSE
)

# Keep only combinations with one variable per category
valid_combinations <- Filter(function(x) check_categories(x, categories_F), combinations)

# ============================================================
# === 4. Compute KMO and Bartlett ===========================
# ============================================================
kmo_scores <- c()
valid_sets <- list()

for (comb in valid_combinations) {
  subset <- na.omit(DataF[, comb, drop=FALSE])
  
  # Skip if too few observations
  if (nrow(subset) < length(comb) + 1) next
  
  # Bartlett test
  bartlett_p <- tryCatch(test_bartlett(subset, comb), error=function(e) NA)
  if (is.na(bartlett_p) || bartlett_p >= 0.05) next
  
  # KMO test
  kmo_val <- tryCatch(calculate_kmo(subset, comb), error=function(e) NA)
  if (!is.na(kmo_val)) {
    kmo_scores <- c(kmo_scores, kmo_val)
    valid_sets <- append(valid_sets, list(comb))
  }
}

# ============================================================
# === 5. Select best combination ============================
# ============================================================
if (length(kmo_scores) > 0) {
  best_index <- which.max(kmo_scores)
  best_combination <- valid_sets[[best_index]]   # ✅ FIX: was best_sets (typo)
  best_kmo <- max(kmo_scores)
  
  cat("✅ Best female variable combination:", paste(best_combination, collapse=", "), "\n")
  cat("✅ Best KMO score:", round(best_kmo, 4), "\n")
} else {
  stop("❌ No valid combination found (Bartlett test not significant or KMO failed).")
}

#The ACP is the same with or without the constraint of reproduction.


# ============================================================
# === 1. Load data and define variables =====================
# ============================================================
Data<-read.csv2("data.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)

# Subset only males
DataM <- subset(Data, Sex == "M")

variables_M <- c(
  "Growth_rate", "Delta_telomere_experiment", "Telomere_pre_hibernation",
  "Temperature_torpor", "Torpor_bout_duration", "Time_inter_torpor", "Time_torpor",
  "log_Mean_cortisol", "Offspring_per_litter", "Offspring_growth_rate",
  "Offspring_number", "Number_transition"
)

categories_M <- list(
  Reproduction = c("Offspring_per_litter", "Offspring_growth_rate", "Offspring_number"),
  Hibernation  = c("Temperature_torpor", "Torpor_bout_duration", "Time_inter_torpor", "Time_torpor"),
  Telomere     = c("Delta_telomere_experiment", "Telomere_pre_hibernation"),
  Cortisol     = c("log_Mean_cortisol"),
  Transition   = c("Number_transition"),
  Growth       = c("Growth_rate")
)

# ============================================================
# === 2. Utility functions ==================================
# ============================================================
calculate_kmo <- function(data, vars) {
  kmo_result <- KMO(data[, vars, drop=FALSE])
  return(kmo_result$MSA[1])
}

test_bartlett <- function(data, vars) {
  corr_mat <- cor(data[, vars, drop=FALSE])
  if (det(corr_mat) < 1e-10) return(1)
  bartlett_test <- cortest.bartlett(corr_mat, n=nrow(data))
  return(bartlett_test$p.value)
}

check_categories <- function(vars, categories) {
  used_cats <- sapply(vars, function(v) {
    cat_name <- names(Filter(function(x) v %in% x, categories))
    if (length(cat_name) == 0) return(NA)
    return(cat_name)
  })
  return(length(unique(used_cats)) == length(vars))
}

# ============================================================
# === 3. Generate all valid combinations ====================
# ============================================================
vars_to_check <- intersect(variables_M, names(DataM))

# Generate all combinations of at least 2 variables
combinations <- unlist(
  lapply(2:length(vars_to_check), function(i) combn(vars_to_check, i, simplify=FALSE)),
  recursive=FALSE
)

# Keep only combinations with one variable per category
valid_combinations <- Filter(function(x) check_categories(x, categories_M), combinations)

# ============================================================
# === 4. Compute KMO and Bartlett ===========================
# ============================================================
kmo_scores <- c()
valid_sets <- list()

for (comb in valid_combinations) {
  subset <- na.omit(DataM[, comb, drop=FALSE])
  
  # Skip if not enough data
  if (nrow(subset) < length(comb) + 1) next
  
  # Bartlett test (must be significant)
  bartlett_p <- tryCatch(test_bartlett(subset, comb), error=function(e) NA)
  if (is.na(bartlett_p) || bartlett_p >= 0.05) next
  
  # KMO test
  kmo_val <- tryCatch(calculate_kmo(subset, comb), error=function(e) NA)
  if (!is.na(kmo_val)) {
    kmo_scores <- c(kmo_scores, kmo_val)
    valid_sets <- append(valid_sets, list(comb))
  }
}

# ============================================================
# === 5. Select best combination ============================
# ============================================================
if (length(kmo_scores) > 0) {
  best_index <- which.max(kmo_scores)
  best_combination <- valid_sets[[best_index]]
  best_kmo <- max(kmo_scores)
  
  cat("✅ Best male variable combination:", paste(best_combination, collapse=", "), "\n")
  cat("✅ Best KMO score:", round(best_kmo, 4), "\n")
} else {
  stop("❌ No valid combination found for males (Bartlett test not significant or KMO failed).")
}


#The ACP is different with and without the constraint of reproduction.

# ============================================================
# === 6. PCA and PC1 extraction =================
# ============================================================
vars_pca <- best_combination
data_pca <- na.omit(DataM[, vars_pca, drop=FALSE])
data_scaled <- scale(data_pca)

res_pca <- PCA(data_scaled, scale.unit=FALSE, ncp=2, graph=FALSE)
DataM$Syndrome_PC1 <- NA
DataM[rownames(data_pca), "Syndrome_PC1"] <- res_pca$ind$coord[,1]
rownames(res_pca$var$coord) <- c("Telomere length variation", "log(Mean cortisol)", "Number of transitions") 

fviz_pca_var(
  res_pca, 
  axes = c(1, 2),
  col.var = "red",
  repel = TRUE,
  arrowsize = 1.2,        
  labelsize = 6,          
  title = "Biplot PCA: variables only",
  circle = TRUE
) + 
  coord_equal() +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_line(color = "grey90"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.margin = margin(10, 10, 10, 10)
  )


# ============================================================
# === 7. PCA robustness check via permutation test ==========
# ============================================================
n_perm <- 500
n_pc <- 2
perm_var <- matrix(NA, nrow=n_perm, ncol=n_pc)
var_obs <- res_pca$eig[1:n_pc,2]

for(p in 1:n_perm) {
  perm_data <- apply(data_scaled, 2, sample)
  res_perm <- tryCatch(PCA(perm_data, scale.unit=FALSE, ncp=n_pc, graph=FALSE), error=function(e) NULL)
  if(!is.null(res_perm)) perm_var[p,] <- res_perm$eig[1:n_pc,2]
}

p_values <- sapply(1:n_pc, function(k) mean(perm_var[,k] >= var_obs[k], na.rm=TRUE))

for(k in 1:n_pc) {
  cat(sprintf("PC%d: Observed variance = %.2f%%, Permutation p-value = %.3f\n", k, var_obs[k], p_values[k]))
}

# ============================================================
# === 8. LME for PC1 ~ Littersize, Birth, Sex =============
# ============================================================
Datana <- na.omit(DataM[, c("Syndrome_PC1","Littersize","Birth","Mother","Sex")])

model_full <- lme(Syndrome_PC1 ~ Littersize + Birth,
                  random = ~1 | Mother, data=Datana, method="ML")

# Model selection using dredge (AICc)
model_selection <- dredge(model_full, fixed = ~ +(1|Mother), rank="AICc", m.lim=c(NA,4))







# ==============================================================================
# =Link between Growth rate and body mass at birth =============================
# ==============================================================================



Datana <- na.omit(Data[, c("Body_mass_birth","Growth_rate","Mother","Sex")])
model <- lme(Growth_rate~Body_mass_birth+Sex,
                  random = ~1 | Mother, data=Datana)
summary(model)
residuals <- resid(model)
plot(density(residuals), main="Residual Density", xlab="Residuals")
curve(dnorm(x, mean=mean(residuals), sd=sd(residuals)), col="red", lwd=2, add=TRUE)
bptest(model_simple)
performance::r2(model_simple)

ggplot(Datana, aes(x=Body_mass_birth, y=Growth_rate, color=Sex)) +
  geom_point(size=3, alpha=0.8) +
  geom_smooth(method="lm", se=TRUE, aes(fill=Sex)) +
  scale_color_manual(values=c("F"="#E41A1C", "M"="#377EB8")) +
  scale_fill_manual(values=c("F"="#E41A1C", "M"="#377EB8")) +
  theme_classic(base_size=14) +
  labs(x="Body mass at birth", y="Growth rate", color="Sex", fill="Sex") +
  theme(legend.position="top")




