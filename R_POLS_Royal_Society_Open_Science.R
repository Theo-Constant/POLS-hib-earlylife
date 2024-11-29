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


#########################
#Test of the POLS hypothesis with all individuals 
##########################
#Variable selection

# Variables
variables <- c("Growth_rate", 
               "Delta_telomere_experiment", "Telomere_pre_hibernation", 
               "Temperature_torpor", 
               "log_Mean_cortisol", "Torpor_bout_duration", 
               "Time_inter_torpor","Time_torpor", "Offspring_per_litter", 
               "Offspring_growth_rate", "Offspring_number","Number_transition")

# Results table
resultats_table <- data.frame(
  Model_Number = integer(),
  N_data = integer(),
  Dependent_Variable = character(),
  Variables = character(),
  Estimate = numeric(),
  Std_Error = numeric(),
  p_value = numeric(),
  Significance = character(),
  Adjusted_p_value = numeric(),
  Adjusted_Significance = character(),
  stringsAsFactors = FALSE
)

# Initialize a global counter for model numbers
model_counter <- 0

# Create a list to store models
model_list <- list()

# Test all variable pairs
for (var1 in variables) {
  for (var2 in variables) {
    # Avoid testing a variable against itself
    if (var1 != var2) {
      
      # Identify necessary columns
      colonnes_utiles <- c(var1, var2, "Mother")
      
      # Filter data to exclude NA in useful columns
      Data_filtered <- na.omit(Data[, colonnes_utiles])
      
      if (nrow(Data_filtered) == 0) {
        next
      }
      
      # Define the model with lme for each pair of variables
      formule <- as.formula(paste(var1, "~", var2))
      model <- lme(formule, random = ~1 | Mother, method = "ML", na.action = "na.fail", data = Data_filtered)
      
      # Increment model counter
      model_counter <- model_counter + 1
      
      # Save the model to the list
      model_list[[paste0("Model_", model_counter)]] <- model
      
      # Extract information on coefficients
      model_summary <- summary(model)
      coef_table <- as.data.frame(model_summary$tTable)
      for (i in 1:nrow(coef_table)) {
        coef_info <- coef_table[i, ]
        n_data <- nrow(Data_filtered)
        
        var_name <- rownames(coef_table)[i]  
        estimate <- as.numeric(coef_info["Value"])
        std_error <- as.numeric(coef_info["Std.Error"])
        Pr_z <- as.numeric(coef_info["p-value"])
        
        # Determine significance level
        if (Pr_z < 0.001) {
          significance <- "***"
        } else if (Pr_z < 0.01) {
          significance <- "**"
        } else if (Pr_z < 0.05) {
          significance <- "*"
        } else {
          significance <- "ns"  
        }
        
        # Add information to the results table
        resultats_table <- rbind(resultats_table, data.frame(
          Model_Number = model_counter,  
          N_data = n_data,
          Dependent_Variable = var1,
          Variables = var_name,
          Estimate = format(round(estimate, 4), nsmall = 3, scientific = FALSE),
          Std_Error = format(round(std_error, 4), nsmall = 3, scientific = FALSE),
          p_value = format(round(Pr_z, 4), nsmall = 3, scientific = FALSE),
          Significance = significance  
        ))
      }
    }
  }
}

# Filter models with at least one significant variable other than the intercept
resultats_significatifs <- resultats_table %>%
  group_by(Model_Number) %>%
  filter(any(Significance %in% c("*", "**", "***") & Variables != "(Intercept)"))

# Generate a PDF file for diagnostics
pdf("Model_Diagnostics.pdf", width = 8, height = 10)

# Apply diagnostics tests to all models
for (model_name in names(model_list)) {
  model <- model_list[[model_name]]
  
  # Add a title page for each model
  plot.new()
  title(main = paste("Diagnostics for", model_name))
  
  # Model summary
  plot.new()
  text(0, 1, paste("Model Summary for", model_name), adj = 0, cex = 1.2)
  model_summary <- capture.output(summary(model))
  text(0, 0.9, paste(model_summary, collapse = "\n"), adj = 0, cex = 0.7)
  
  # Density Plot of Residuals with Normal Curve Overlay
  residuals <- resid(model)
  plot(density(residuals), main = paste("Density Plot of Residuals for", model_name), xlab = "Residuals", col = "blue", lwd = 2)
  curve(dnorm(x, mean = mean(residuals), sd = sd(residuals)), col = "red", lwd = 2, add = TRUE)
  legend("topright", legend = c("Residuals Density", "Normal Distribution"), col = c("blue", "red"), lwd = c(2, 2))
  
  # QQ-Plot of Residuals
  qqPlot(residuals, main = paste("QQ-Plot of Residuals for", model_name))  # Requires car::qqPlot
  
  # Breusch-Pagan test
  bptest_result <- bptest(model)  # Requires lmtest::bptest
  plot.new()
  text(0, 1, paste("Breusch-Pagan Test Result for", model_name), adj = 0, cex = 1.2)
  bptest_text <- capture.output(bptest_result)
  text(0, 0.9, paste(bptest_text, collapse = "\n"), adj = 0, cex = 0.7)
}

# Close the PDF device
dev.off()

# Show significant results
View(resultats_table)
View(resultats_significatifs)

###############
# Create a Word document for supplementary materials
doc <- read_docx()

# Add a title for the complete table
doc <- body_add_par(doc, value = "Results Table (All Models)", style = "heading 1")

# Create a flextable for results_table
flextable_resultats <- flextable(resultats_table) %>%
  theme_vanilla() %>%
  autofit()  

# Add complete table to Word document
doc <- body_add_flextable(doc, flextable_resultats)

# Add a page break
doc <- body_add_par(doc, "", style = "Normal")
doc <- body_add_par(doc, " ", style = "Normal")
doc <- body_add_par(doc, "Results Table (Significant Models)", style = "heading 1")

# Create a flextable for results_significant
flextable_significatifs <- flextable(resultats_significatifs) %>%
  theme_vanilla() %>%
  autofit()  

# Add the table of significant results to the Word document
doc <- body_add_flextable(doc, flextable_significatifs)

# Save Word document
output_file <- "Results_Tables.docx"
print(doc, target = output_file)


###########
#new selection of variable for test ACP based on KMO value

dataACP1 <- Data[,c("Offspring_number","Offspring_growth_rate","Number_transition","log_Mean_cortisol","Delta_telomere_experiment","Offspring_per_litter","Time_torpor")]
imputed_data1 <- imputePCA(dataACP1)
dataimp1<-imputed_data1$completeObs
dataimp1<-as.data.frame(dataimp1)

# Function to calculate the KMO for a subset of variables
calculate_kmo <- function(data, vars) {
  kmo_result <- KMO(data[, vars])
  return(kmo_result$MSA)  
}

# Bartlett test function
test_bartlett <- function(data, vars) {
  correlation_matrix <- cor(data[, vars])
  bartlett_test <- cortest.bartlett(correlation_matrix, n = nrow(data))
  return(bartlett_test$p.value) 
}

# Test all possible variable combinations (from 3 to 6 variables)
combinations <- unlist(lapply(3:ncol(dataimp1), function(i) combn(names(dataimp1), i, simplify = FALSE)), recursive = FALSE)

# Initialize lists to store results
kmo_scores <- c()
valid_combinations <- list()

# Calculate KMO only for combinations with a significant Bartlett test
for (comb in combinations) {
  bartlett_pvalue <- test_bartlett(dataimp1, comb)
  
  if (bartlett_pvalue < 0.05) {  
    kmo_value <- calculate_kmo(dataimp1, comb)  
    kmo_scores <- c(kmo_scores, kmo_value)
    valid_combinations <- append(valid_combinations, list(comb))
  }
}

# Find the combination that gives the best KMO among the valid ones
if (length(kmo_scores) > 0) {
  best_combination <- valid_combinations[[which.max(kmo_scores)]]
  best_kmo <- max(kmo_scores)
  
  # Display result 
  cat("The best combination of variables is :", paste(best_combination, collapse = ", "), "\n")
  cat("The best KMO score (Overall MSA) is :", best_kmo, "\n")
} else {
  cat("No significant combination was found with a Bartlett test.")
}

##################
#ACP on selected variable

dataACP2 <- Data[,c("Offspring_growth_rate","Number_transition","Delta_telomere_experiment","Offspring_per_litter","Time_torpor")]
imputed_data2 <- imputePCA(dataACP2)
dataimp2<-imputed_data2$completeObs
nrow(dataimp2)
R2<-cor(dataimp2)
cortest.bartlett(R2,n=34)
KMO(R2)
dataimp2<-as.data.frame(dataimp2)
is.data.frame(dataimp2)

#############
#permutation test

permutation_test_acp <- function(data, n_permutations = 1000) {
  if (!is.data.frame(data)) stop("The data must be a data.frame.")
  if (anyNA(data)) stop("The dataset contains missing values.")
  
  # Perform PCA on original data
  res_pca <- PCA(data, scale.unit = TRUE, graph = FALSE)
  obs_values <- res_pca$eig[, 1]  
  num_pcs <- length(obs_values)
  
  # Initialization to store permutation results
  permuted_values <- matrix(NA, ncol = num_pcs, nrow = n_permutations)
  
  for (i in 1:n_permutations) {
    permuted_data <- as.data.frame(lapply(data, sample))
    
    perm_res <- tryCatch(
      PCA(permuted_data, scale.unit = TRUE, graph = FALSE),
      error = function(e) {
        message(paste("Error during PCA at iteration", i, ":", e$message))
        return(NULL)
      }
    )
    
    if (is.null(perm_res) || nrow(perm_res$eig) != num_pcs) next
    
    permuted_values[i, ] <- perm_res$eig[, 1]
  }
  
  print("Dimensions of swapped values :")
  print(dim(permuted_values))
  print("Length of observed eigenvalues :")
  print(length(obs_values))
  
  # Calculating p-values
  p_values <- sapply(1:ncol(permuted_values), function(j) {
    mean(permuted_values[, j] >= obs_values[j], na.rm = TRUE)
  })
  
  results <- data.frame(
    PC = paste0("PC", 1:num_pcs),
    Observed_Eigenvalue = obs_values,
    p_value = p_values
  )
  
  return(results)
}

# Apply the permutation test function to a dataset
res_perm <- permutation_test_acp(dataimp2, n_permutations = 1000)
print(res_perm)

##################
#evalution of ACP 

res.PCA<-PCA(dataimp2, graph = TRUE)
eig.val <- get_eigenvalue(res.PCA)
eig.val
# Extract variable contributions for each axis
contributions <- res.PCA$var$contrib
print(contributions)

#extract the first component and test effect of sex, littersize and Birth
Data$ACPall<-res.PCA$ind$coord[,1]
test1<-lme(ACPall ~Littersize*Sex+Birth*Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(model.avg(dredge(test1, fixed=~+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test1<-lme(ACPall ~Littersize*Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(test1)

#validation
residuals <- resid(test1)
#density
plot(density(residuals), main = paste("Density Plot of Residuals for", model_name), xlab = "Residuals", col = "blue", lwd = 2)
curve(dnorm(x, mean = mean(residuals), sd = sd(residuals)), col = "red", lwd = 2, add = TRUE)
legend("topright", legend = c("Residuals Density", "Normal Distribution"), col = c("blue", "red"), lwd = c(2, 2))

# QQ-Plot of Residuals
qqPlot(residuals, main = paste("QQ-Plot of Residuals for", model_name))  # Requires car::qqPlot

# Breusch-Pagan test
bptest(test1)

# Calculate R?
r2_results <- r2(test1)
print(r2_results)

#plot
ggplot(Data, aes(x = ACPall, y = Littersize, color = Sex)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", aes(group = Sex, color = Sex), se = FALSE, linetype = "solid") +
  theme_classic() +
  labs(x = "PC1", y = "Litter size")

######################
###Test only on female

# Variables
variables <- c("Growth_rate", 
               "Delta_telomere_experiment", "Telomere_pre_hibernation", 
               "Temperature_torpor", 
               "log_Mean_cortisol", "Torpor_bout_duration", 
               "Time_inter_torpor","Time_torpor", "Offspring_per_litter", 
               "Offspring_growth_rate", "Offspring_number","Number_transition")

# Results table
resultats_table <- data.frame(
  Model_Number = integer(),
  N_data = integer(),
  Dependent_Variable = character(),
  Variables = character(),
  Estimate = numeric(),
  Std_Error = numeric(),
  p_value = numeric(),
  Significance = character(),
  Adjusted_p_value = numeric(),
  Adjusted_Significance = character(),
  stringsAsFactors = FALSE
)

# Initialize a global counter for model numbers
model_counter <- 0

# Create a list to store models
model_list <- list()

# Test all variable pairs
for (var1 in variables) {
  for (var2 in variables) {
    # Avoid testing a variable against itself
    if (var1 != var2) {
      
      # Identify necessary columns
      colonnes_utiles <- c(var1, var2, "Mother")
      
      # Filter data to exclude NA in useful columns
      DataF<-subset(Data,Data$Sex=="F")
      Data_filtered <- na.omit(DataF[, colonnes_utiles])
      
      if (nrow(Data_filtered) == 0) {
        next
      }
      
      # Define the model with lme for each pair of variables
      formule <- as.formula(paste(var1, "~", var2))
      model <- lme(formule, random = ~1 | Mother, method = "ML", na.action = "na.fail", data = Data_filtered)
      
      # Increment model counter
      model_counter <- model_counter + 1
      
      # Save the model to the list
      model_list[[paste0("Model_", model_counter)]] <- model
      
      # Extract information on coefficients
      model_summary <- summary(model)
      coef_table <- as.data.frame(model_summary$tTable)
      for (i in 1:nrow(coef_table)) {
        coef_info <- coef_table[i, ]
        n_data <- nrow(Data_filtered)
        
        var_name <- rownames(coef_table)[i]  
        estimate <- as.numeric(coef_info["Value"])
        std_error <- as.numeric(coef_info["Std.Error"])
        Pr_z <- as.numeric(coef_info["p-value"])
        
        # Determine significance level
        if (Pr_z < 0.001) {
          significance <- "***"
        } else if (Pr_z < 0.01) {
          significance <- "**"
        } else if (Pr_z < 0.05) {
          significance <- "*"
        } else {
          significance <- "ns"  
        }
        
        # Add information to the results table
        resultats_table <- rbind(resultats_table, data.frame(
          Model_Number = model_counter,  
          N_data = n_data,
          Dependent_Variable = var1,
          Variables = var_name,
          Estimate = format(round(estimate, 4), nsmall = 3, scientific = FALSE),
          Std_Error = format(round(std_error, 4), nsmall = 3, scientific = FALSE),
          p_value = format(round(Pr_z, 4), nsmall = 3, scientific = FALSE),
          Significance = significance  
        ))
      }
    }
  }
}

# Filter models with at least one significant variable other than the intercept
resultats_significatifs <- resultats_table %>%
  group_by(Model_Number) %>%
  filter(any(Significance %in% c("*", "**", "***") & Variables != "(Intercept)"))

# Generate a PDF file for diagnostics
pdf("Model_Diagnostics.pdf", width = 8, height = 10)

# Apply diagnostics tests to all models
for (model_name in names(model_list)) {
  model <- model_list[[model_name]]
  
  # Add a title page for each model
  plot.new()
  title(main = paste("Diagnostics for", model_name))
  
  # Model summary
  plot.new()
  text(0, 1, paste("Model Summary for", model_name), adj = 0, cex = 1.2)
  model_summary <- capture.output(summary(model))
  text(0, 0.9, paste(model_summary, collapse = "\n"), adj = 0, cex = 0.7)
  
  # Density Plot of Residuals with Normal Curve Overlay
  residuals <- resid(model)
  plot(density(residuals), main = paste("Density Plot of Residuals for", model_name), xlab = "Residuals", col = "blue", lwd = 2)
  curve(dnorm(x, mean = mean(residuals), sd = sd(residuals)), col = "red", lwd = 2, add = TRUE)
  legend("topright", legend = c("Residuals Density", "Normal Distribution"), col = c("blue", "red"), lwd = c(2, 2))
  
  # QQ-Plot of Residuals
  qqPlot(residuals, main = paste("QQ-Plot of Residuals for", model_name))  # Requires car::qqPlot
  
  # Breusch-Pagan test
  bptest_result <- bptest(model)  # Requires lmtest::bptest
  plot.new()
  text(0, 1, paste("Breusch-Pagan Test Result for", model_name), adj = 0, cex = 1.2)
  bptest_text <- capture.output(bptest_result)
  text(0, 0.9, paste(bptest_text, collapse = "\n"), adj = 0, cex = 0.7)
}

# Close the PDF device
dev.off()

# Show significant results
View(resultats_table)
View(resultats_significatifs)

###############
# Create a Word document for supplementary materials
doc <- read_docx()

# Add a title for the complete table
doc <- body_add_par(doc, value = "Results Table (All Models)", style = "heading 1")

# Create a flextable for results_table
flextable_resultats <- flextable(resultats_table) %>%
  theme_vanilla() %>%
  autofit()  

# Add complete table to Word document
doc <- body_add_flextable(doc, flextable_resultats)

# Add a page break
doc <- body_add_par(doc, "", style = "Normal")
doc <- body_add_par(doc, " ", style = "Normal")
doc <- body_add_par(doc, "Results Table (Significant Models)", style = "heading 1")

# Create a flextable for results_significant
flextable_significatifs <- flextable(resultats_significatifs) %>%
  theme_vanilla() %>%
  autofit()  

# Add the table of significant results to the Word document
doc <- body_add_flextable(doc, flextable_significatifs)

# Save Word document
output_file <- "Results_Tables.docx"
print(doc, target = output_file)

##########################
#new variable selection based on KMO 
dataACPF1 <- DataF[,c("log_Mean_cortisol","Delta_telomere_experiment","Offspring_per_litter","Number_transition","Growth_rate","Offspring_number","Time_inter_torpor","Offspring_growth_rate")]
imputed_dataF1 <- imputePCA(dataACPF1)
dataimpF1<-imputed_dataF1$completeObs
dataimpF1<-as.data.frame(dataimpF1)
dataimpF1

# Function to calculate the KMO for a subset of variables
calculate_kmo <- function(data, vars) {
  kmo_result <- KMO(data[, vars])
  return(kmo_result$MSA)  
}

# Bartlett test function
test_bartlett <- function(data, vars) {
  correlation_matrix <- cor(data[, vars])
  bartlett_test <- cortest.bartlett(correlation_matrix, n = nrow(data))
  return(bartlett_test$p.value) 
}

# Test all possible variable combinations (from 3 to 6 variables)
combinations <- unlist(lapply(3:ncol(dataimpF1), function(i) combn(names(dataimpF1), i, simplify = FALSE)), recursive = FALSE)

# Initialize lists to store results
kmo_scores <- c()
valid_combinations <- list()

# Calculate KMO only for combinations with a significant Bartlett test
for (comb in combinations) {
  bartlett_pvalue <- test_bartlett(dataimpF1, comb)
  
  if (bartlett_pvalue < 0.05) {  
    kmo_value <- calculate_kmo(dataimpF1, comb)  
    kmo_scores <- c(kmo_scores, kmo_value)
    valid_combinations <- append(valid_combinations, list(comb))
  }
}

# Find the combination that gives the best KMO among the valid ones
if (length(kmo_scores) > 0) {
  best_combination <- valid_combinations[[which.max(kmo_scores)]]
  best_kmo <- max(kmo_scores)
  
  # Display result 
  cat("The best combination of variables is :", paste(best_combination, collapse = ", "), "\n")
  cat("The best KMO score (Overall MSA) is :", best_kmo, "\n")
} else {
  cat("No significant combination was found with a Bartlett test.")
}

#ACP test on selected variable

dataACPF2 <- DataF[,c("Offspring_per_litter","Number_transition","Growth_rate","Offspring_growth_rate")]
imputed_dataF2 <- imputePCA(dataACPF2)
dataimpF2<-imputed_dataF2$completeObs
nrow(dataimpF2)
R2<-cor(dataimpF2)
cortest.bartlett(R2,n=16)
KMO(R2)

#Evoluation of PCA
res.PCAF1<-PCA(dataimpF2, graph = TRUE)
eig.val <- get_eigenvalue(res.PCAF1)
eig.val

# Extraire les contributions des variables pour chaque axe
contributions <- res.PCA1$var$contrib
print(contributions)

#permutation test
permutation_test_acp <- function(data, n_permutations = 1000) {
  if (!is.data.frame(data)) stop("Les données doivent être un data.frame.")
  if (anyNA(data)) stop("Le jeu de données contient des valeurs manquantes.")
  
  # Réaliser l'ACP sur les données originales
  res_pca <- PCA(data, scale.unit = TRUE, graph = FALSE)
  obs_values <- res_pca$eig[, 1]  # Valeurs propres observées
  num_pcs <- length(obs_values)
  
  # Initialisation pour stocker les résultats des permutations
  permuted_values <- matrix(NA, ncol = num_pcs, nrow = n_permutations)
  
  for (i in 1:n_permutations) {
    permuted_data <- as.data.frame(lapply(data, sample))
    
    perm_res <- tryCatch(
      PCA(permuted_data, scale.unit = TRUE, graph = FALSE),
      error = function(e) {
        message(paste("Erreur lors de l'ACP à l'itération", i, ":", e$message))
        return(NULL)
      }
    )
    
    if (is.null(perm_res) || nrow(perm_res$eig) != num_pcs) next
    
    permuted_values[i, ] <- perm_res$eig[, 1]
  }
  
  print("Dimensions des valeurs permutées :")
  print(dim(permuted_values))
  print("Longueur des valeurs propres observées :")
  print(length(obs_values))
  
  # Calcul des p-valeurs
  p_values <- sapply(1:ncol(permuted_values), function(j) {
    mean(permuted_values[, j] >= obs_values[j], na.rm = TRUE)
  })
  
  results <- data.frame(
    PC = paste0("PC", 1:num_pcs),
    Observed_Eigenvalue = obs_values,
    p_value = p_values
  )
  
  return(results)
}


# Apply the permutation test function to a dataset
dataimpF2<-as.data.frame(dataimpF2)
is.data.frame(dataimpF2)
res_permF <- permutation_test_acp(dataimpF2, n_permutations = 1000)
print(res_permF)

#test the effect of littersize, sex and birth on principal component
DataF$ACPF<-res.PCAF1$ind$coord[,1]
test1F<-lme(ACPF ~Littersize+Birth, random=~1|Mother,method="ML", na.action = "na.fail", data=DataF)
summary(model.avg(dredge(test1F, fixed=~+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test1F<-lme(ACPF ~Littersize, random=~1|Mother,method="ML", na.action = "na.fail", data=DataF)
summary(test1F)

#validation
residuals <- resid(test1F)
#density
plot(density(residuals), main = paste("Density Plot of Residuals for", model_name), xlab = "Residuals", col = "blue", lwd = 2)
curve(dnorm(x, mean = mean(residuals), sd = sd(residuals)), col = "red", lwd = 2, add = TRUE)
legend("topright", legend = c("Residuals Density", "Normal Distribution"), col = c("blue", "red"), lwd = c(2, 2))

# QQ-Plot of Residuals
qqPlot(residuals, main = paste("QQ-Plot of Residuals for", model_name))  # Requires car::qqPlot

# Breusch-Pagan test
bptest(test1F)

# Calculate R?
r2_results <- r2(test1F)
print(r2_results)

#plot
ggplot(DataF, aes(x = ACPF, y = Littersize)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = "PC1F", y = "Litter size")


##################
###Test only on male


variables <- c("Growth_rate", 
               "Delta_telomere_experiment", "Telomere_pre_hibernation", 
               "Temperature_torpor", 
               "log_Mean_cortisol", "Torpor_bout_duration", 
               "Time_inter_torpor","Time_torpor", "Offspring_per_litter", 
               "Offspring_growth_rate", "Offspring_number","Number_transition","Testosterone_post_hibernation")

# Results table
resultats_table <- data.frame(
  Model_Number = integer(),
  N_data = integer(),
  Dependent_Variable = character(),
  Variables = character(),
  Estimate = numeric(),
  Std_Error = numeric(),
  p_value = numeric(),
  Significance = character(),
  Adjusted_p_value = numeric(),
  Adjusted_Significance = character(),
  stringsAsFactors = FALSE
)

# Initialize a global counter for model numbers
model_counter <- 0

# Create a list to store models
model_list <- list()

# Test all variable pairs
for (var1 in variables) {
  for (var2 in variables) {
    # Avoid testing a variable against itself
    if (var1 != var2) {
      
      # Identify necessary columns
      colonnes_utiles <- c(var1, var2, "Mother")
      
      # Filter data to exclude NA in useful columns
      DataM<-subset(Data,Data$Sex=="M")
      Data_filtered <- na.omit(DataM[, colonnes_utiles])
      
      if (nrow(Data_filtered) == 0) {
        next
      }
      
      # Define the model with lme for each pair of variables
      formule <- as.formula(paste(var1, "~", var2))
      model <- lme(formule, random = ~1 | Mother, method = "ML", na.action = "na.fail", data = Data_filtered)
      
      # Increment model counter
      model_counter <- model_counter + 1
      
      # Save the model to the list
      model_list[[paste0("Model_", model_counter)]] <- model
      
      # Extract information on coefficients
      model_summary <- summary(model)
      coef_table <- as.data.frame(model_summary$tTable)
      for (i in 1:nrow(coef_table)) {
        coef_info <- coef_table[i, ]
        n_data <- nrow(Data_filtered)
        
        var_name <- rownames(coef_table)[i]  
        estimate <- as.numeric(coef_info["Value"])
        std_error <- as.numeric(coef_info["Std.Error"])
        Pr_z <- as.numeric(coef_info["p-value"])
        
        # Determine significance level
        if (Pr_z < 0.001) {
          significance <- "***"
        } else if (Pr_z < 0.01) {
          significance <- "**"
        } else if (Pr_z < 0.05) {
          significance <- "*"
        } else {
          significance <- "ns"  
        }
        
        # Add information to the results table
        resultats_table <- rbind(resultats_table, data.frame(
          Model_Number = model_counter,  
          N_data = n_data,
          Dependent_Variable = var1,
          Variables = var_name,
          Estimate = format(round(estimate, 4), nsmall = 3, scientific = FALSE),
          Std_Error = format(round(std_error, 4), nsmall = 3, scientific = FALSE),
          p_value = format(round(Pr_z, 4), nsmall = 3, scientific = FALSE),
          Significance = significance  
        ))
      }
    }
  }
}

# Filter models with at least one significant variable other than the intercept
resultats_significatifs <- resultats_table %>%
  group_by(Model_Number) %>%
  filter(any(Significance %in% c("*", "**", "***") & Variables != "(Intercept)"))

# Generate a PDF file for diagnostics
pdf("Model_Diagnostics.pdf", width = 8, height = 10)

# Apply diagnostics tests to all models
for (model_name in names(model_list)) {
  model <- model_list[[model_name]]
  
  # Add a title page for each model
  plot.new()
  title(main = paste("Diagnostics for", model_name))
  
  # Model summary
  plot.new()
  text(0, 1, paste("Model Summary for", model_name), adj = 0, cex = 1.2)
  model_summary <- capture.output(summary(model))
  text(0, 0.9, paste(model_summary, collapse = "\n"), adj = 0, cex = 0.7)
  
  # Density Plot of Residuals with Normal Curve Overlay
  residuals <- resid(model)
  plot(density(residuals), main = paste("Density Plot of Residuals for", model_name), xlab = "Residuals", col = "blue", lwd = 2)
  curve(dnorm(x, mean = mean(residuals), sd = sd(residuals)), col = "red", lwd = 2, add = TRUE)
  legend("topright", legend = c("Residuals Density", "Normal Distribution"), col = c("blue", "red"), lwd = c(2, 2))
  
  # QQ-Plot of Residuals
  qqPlot(residuals, main = paste("QQ-Plot of Residuals for", model_name))  # Requires car::qqPlot
  
  # Breusch-Pagan test
  bptest_result <- bptest(model)  # Requires lmtest::bptest
  plot.new()
  text(0, 1, paste("Breusch-Pagan Test Result for", model_name), adj = 0, cex = 1.2)
  bptest_text <- capture.output(bptest_result)
  text(0, 0.9, paste(bptest_text, collapse = "\n"), adj = 0, cex = 0.7)
}

# Close the PDF device
dev.off()

# Show significant results
View(resultats_table)
View(resultats_significatifs)

###############
# Create a Word document for supplementary materials
doc <- read_docx()

# Add a title for the complete table
doc <- body_add_par(doc, value = "Results Table (All Models)", style = "heading 1")

# Create a flextable for results_table
flextable_resultats <- flextable(resultats_table) %>%
  theme_vanilla() %>%
  autofit()  

# Add complete table to Word document
doc <- body_add_flextable(doc, flextable_resultats)

# Add a page break
doc <- body_add_par(doc, "", style = "Normal")
doc <- body_add_par(doc, " ", style = "Normal")
doc <- body_add_par(doc, "Results Table (Significant Models)", style = "heading 1")

# Create a flextable for results_significant
flextable_significatifs <- flextable(resultats_significatifs) %>%
  theme_vanilla() %>%
  autofit()  

# Add the table of significant results to the Word document
doc <- body_add_flextable(doc, flextable_significatifs)

# Save Word document
output_file <- "Results_Tables.docx"
print(doc, target = output_file)

##########
#new variable selection based on KMO

dataACPM1 <- DataM[,c("Number_transition","Delta_telomere_experiment","Growth_rate","Telomere_pre_hibernation")]
imputed_dataM1 <- imputePCA(dataACPM1)
dataimpM1<-imputed_dataM1$completeObs
dataimpM1<-as.data.frame(dataACPM1)
dataimpM1

# Function to calculate the KMO for a subset of variables
calculate_kmo <- function(data, vars) {
  kmo_result <- KMO(data[, vars])
  return(kmo_result$MSA)  
}

# Bartlett test function
test_bartlett <- function(data, vars) {
  correlation_matrix <- cor(data[, vars])
  bartlett_test <- cortest.bartlett(correlation_matrix, n = nrow(data))
  return(bartlett_test$p.value) 
}

# Test all possible variable combinations (from 3 to 6 variables)
combinations <- unlist(lapply(3:ncol(dataimpM1), function(i) combn(names(dataimpM1), i, simplify = FALSE)), recursive = FALSE)

# Initialize lists to store results
kmo_scores <- c()
valid_combinations <- list()

# Calculate KMO only for combinations with a significant Bartlett test
for (comb in combinations) {
  bartlett_pvalue <- test_bartlett(dataimpM1, comb)
  
  if (bartlett_pvalue < 0.05) {  
    kmo_value <- calculate_kmo(dataimpM1, comb)  
    kmo_scores <- c(kmo_scores, kmo_value)
    valid_combinations <- append(valid_combinations, list(comb))
  }
}

# Find the combination that gives the best KMO among the valid ones
if (length(kmo_scores) > 0) {
  best_combination <- valid_combinations[[which.max(kmo_scores)]]
  best_kmo <- max(kmo_scores)
  
  # Display result 
  cat("The best combination of variables is :", paste(best_combination, collapse = ", "), "\n")
  cat("The best KMO score (Overall MSA) is :", best_kmo, "\n")
} else {
  cat("No significant combination was found with a Bartlett test.")
}

#ACP test on selected variable
R2<-cor(dataACPM1)
cortest.bartlett(R2,n=16)
KMO(R2)
res.PCAM1<-PCA(dataACPM1, graph = TRUE)
eig.val <- get_eigenvalue(res.PCAM1)
eig.val

# Extraire les contributions des variables pour chaque axe
contributions <- res.PCAM1$var$contrib
print(contributions)

#permutation test
permutation_test_acp <- function(data, n_permutations = 1000) {
  if (!is.data.frame(data)) stop("Les données doivent être un data.frame.")
  if (anyNA(data)) stop("Le jeu de données contient des valeurs manquantes.")
  
  # Réaliser l'ACP sur les données originales
  res_pca <- PCA(data, scale.unit = TRUE, graph = FALSE)
  obs_values <- res_pca$eig[, 1]  # Valeurs propres observées
  num_pcs <- length(obs_values)
  
  # Initialisation pour stocker les résultats des permutations
  permuted_values <- matrix(NA, ncol = num_pcs, nrow = n_permutations)
  
  for (i in 1:n_permutations) {
    permuted_data <- as.data.frame(lapply(data, sample))
    
    perm_res <- tryCatch(
      PCA(permuted_data, scale.unit = TRUE, graph = FALSE),
      error = function(e) {
        message(paste("Erreur lors de l'ACP à l'itération", i, ":", e$message))
        return(NULL)
      }
    )
    
    if (is.null(perm_res) || nrow(perm_res$eig) != num_pcs) next
    
    permuted_values[i, ] <- perm_res$eig[, 1]
  }
  
  print("Dimensions des valeurs permutées :")
  print(dim(permuted_values))
  print("Longueur des valeurs propres observées :")
  print(length(obs_values))
  
  # Calcul des p-valeurs
  p_values <- sapply(1:ncol(permuted_values), function(j) {
    mean(permuted_values[, j] >= obs_values[j], na.rm = TRUE)
  })
  
  results <- data.frame(
    PC = paste0("PC", 1:num_pcs),
    Observed_Eigenvalue = obs_values,
    p_value = p_values
  )
  
  return(results)
}


# Apply the permutation test function to a dataset
dataimpM1<-as.data.frame(dataimpM1)
is.data.frame(dataimpM1)
res_permM <- permutation_test_acp(dataimpM1, n_permutations = 1000)
print(res_permM)

#test the effect of sex, littersize and Birth on principal component
DataM$ACPM<-res.PCAM1$ind$coord[,1]
testM1<-lme(ACPM ~Littersize+Birth, random=~1|Mother,method="ML", na.action = "na.fail", data=DataM)
summary(model.avg(dredge(testM1, fixed=~+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
testM1<-lme(ACPM ~Littersize, random=~1|Mother,method="ML", na.action = "na.fail", data=DataM)
summary(testM1)

#validation
residuals <- resid(testM1)
#density
plot(density(residuals), main = paste("Density Plot of Residuals for", model_name), xlab = "Residuals", col = "blue", lwd = 2)
curve(dnorm(x, mean = mean(residuals), sd = sd(residuals)), col = "red", lwd = 2, add = TRUE)
legend("topright", legend = c("Residuals Density", "Normal Distribution"), col = c("blue", "red"), lwd = c(2, 2))

# QQ-Plot of Residuals
qqPlot(residuals, main = paste("QQ-Plot of Residuals for", model_name))  # Requires car::qqPlot

# Breusch-Pagan test
bptest(testM1)

# Calculate R?
r2_results <- r2(testM1)
print(r2_results)


