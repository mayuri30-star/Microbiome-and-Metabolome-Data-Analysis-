#installing the necessary libraries
install.packages(c("ggplot2","dplyr","readr"))
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

#Loading of the datsets

#Microbiome data
microbiome_data <- read.csv("microbiome_data.csv")

#Metabolome_data
metabolome_data <- read.csv("metabolome_data.csv")

#Metadata
metadata <- read.csv("metadata.csv")

#Displaying the first few rows of each dataset 
head(microbiome_data)
head(metabolome_data)
head(metadata)

#Checking missing values in microbiome,metabolome and metadata
sum(is.na(microbiome_data))
sum(is.na(metabolome_data))
sum(is.na(metadata))

#Merging the microbiome data with metadata into one combined dataset
combine_data <- metadata %>%
  left_join(microbiome_data,by ="SampleID") %>%
  left_join(metabolome_data,by ="SampleID")

#Viewing the first few rows 
head(combine_data)
dim(combine_data)         #cheking the dimentions of merged dataset


#Part 1: Descriptive statistics and data visualization

#Step 1: Descriptive Statistics

#Microbial species summary statistics
microbiome_summary <- microbiome_data %>%
  summarise(across(everything(), list(mean=mean,median=median, sd = sd), na.rm = TRUE))

#Metabolite concentration summary statistics
metabolome_summary <- metabolome_data %>%
  summarise(across(everything(), list(mean=mean, median=median,sd=sd),na.rm=TRUE))

#The Summary statistics
print(microbiome_summary)
print(metabolome_summary)

#Step 2: Data Visualization - Boxplots

# Coverting the 'Group' column to a factor
combine_data$Group <- factor(combine_data$Group, levels = c("Healthy", "Disease"))

# Extracting the names of microbial species removing sample ID
species_names <- colnames(combine_data)[-1]

# Extracting the names of metabolites by removing SampleID
metabolite_names <- colnames(metabolome_data)[-1]


# Boxplot for metabolite concentrations by group

# Creating an empty list for storing metabolite plots
plot_list <- list()   

# Loop through each metabolite and generate boxplots
for (metabolite in metabolite_names) {
  plot <- ggplot(combine_data, aes(x = Group, y = .data[[metabolite]], fill = Group)) +
    geom_boxplot() +
    labs(
      title = paste("Concentration of", metabolite, "by Group"),
      x = "Group",
      y = "Concentration"
    ) +
    theme_minimal()
  
  # Adding each plot to the list
  
  plot_list[[metabolite]] <- plot
}

# Combining all the plots into single layout using patchwork
combined_plot <- wrap_plots(plot_list)+ 
  plot_layout(guides= "collect") &
  theme(legend.position = "bottom")
print(combined_plot)

#Boxplots for microbial species abundance by group

#Creating empty list to store the species plots
plot_list <- list()

# Loop through each species and generate a boxplot
for (species in species_names) {
  plot <- ggplot(combine_data, aes_string(x = "Group", y = species, fill = "Group")) +
    geom_boxplot() +
    labs(
      title = paste("Abundance of", species, "by Group"),
      x = "Group",
      y = "Abundance"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  plot_list[[species]] <- plot
}

# Combining all plots into a single layout
combined_plot <- wrap_plots(plot_list)+ 
  plot_layout(guides= "collect")
print(combined_plot)

# Part 2: T-test and Mann-Whitney-Wilcoxon Test

# Initialize an empty data frame to store results for both tests 
result_df <- data.frame(
  Species = character(),
  Test = character(),
  p_value = numeric(),
  Significance = character(),
  stringsAsFactors = FALSE
)

# Extracting the names of microbial species by removing names
species_names <- colnames(combine_data)[!(colnames(combine_data) %in% c("SampleID", "Group"))]

# Loop through each microbial species and perform T-test and Mann-Whitney-Wilcoxon test
for (species in species_names) {
  
  # Filtering the data that include only valid groups
  clean_data <- combine_data %>%
    filter(Group %in% c("Healthy", "Disease") & !is.na(.data[[species]]))
  
  # Checking if 'Group' has exactly two levels 
  if (length(unique(clean_data$Group)) == 2) {
    
  # Performing an independent t-test for each species
    t_test_result <- t.test(clean_data[[species]] ~ clean_data$Group)
    t_p_value <- t_test_result$p.value
    t_significance <- ifelse(t_p_value < 0.05, "Significant", "Not Significant")
  
  # Append T-test results   
    result_df <- rbind(result_df, data.frame(Species = species, Test = "T-test", p_value = t_p_value, Significance = t_significance))
    
  # Performing Mann-Whitney-Wilcoxon test as a non-parametric alternative
    wilcoxon_test_result <- wilcox.test(clean_data[[species]] ~ clean_data$Group)
    wilcoxon_p_value <- wilcoxon_test_result$p.value
    wilcoxon_significance <- ifelse(wilcoxon_p_value < 0.05, "Significant", "Not Significant")
    
  # Append Mann-Whitney-Wilcoxon test results 
    result_df <- rbind(result_df, data.frame(Species = species, Test = "Mann-Whitney-Wilcoxon", p_value = wilcoxon_p_value, Significance = wilcoxon_significance))
    
  } else {
    
  # Skip the test if 'Group' does not have exactly two levels
    cat("\nSkipping", species, "because 'Group' does not have exactly two levels.\n")
  }
  
  cat("\n---------------------\n")          # The Seperator for readability
}

# Display the result for both tests
print(result_df)


# Part 3: ANOVA and Kruskal-Wallis Tests on Metabollites

# Loading the necessary libraries
library(dplyr)

# Displaying the first few rows of the loaded dataset
cat("First few lines of the full dataset:\n")
print(head(metabolome_data))

# Merging with metadata to include the group information
metabolome_data <- metabolome_data %>%
  left_join(metadata %>% select(SampleID, Group), by = "SampleID")

# Selecting only 'Group' and metabolites columns by removing rows with missing values
metabolome_data <- metabolome_data %>%
  select(Group, starts_with("Metabolite")) %>%
  na.omit()

# An empty dataframe to store the test results 
result_df <- data.frame(
  Metabolite = character(),
  Test = character(),
  p_value = numeric(),
  Significance = character(),
  stringsAsFactors = FALSE
)

#  Make a loop through each metabolite ( from Metabolite_1 to Metabolite_10)
for (metabolite in paste0("Metabolite_", 1:10)) {
  
  # Extracting the data
  metabolite_data <- metabolome_data %>%
    select(Group, all_of(metabolite))
  
# Performing the ANOVA test to check if there is a significant diffrence between groups
  anova_result <- aov(as.formula(paste(metabolite, "~ Group")), data = metabolite_data)
  anova_p_value <- summary(anova_result)[[1]]$`Pr(>F)`[1] 
  anova_significance <- ifelse(anova_p_value < 0.05, "Significant", "Not Significant")
  
# Append ANOVA results to the result data frame
  result_df <- rbind(result_df, data.frame(
    Metabolite = metabolite, 
    Test = "ANOVA", 
    p_value = anova_p_value, 
    Significance = anova_significance
  ))
  
# Performing Kruskal-Wallis test (non-parametric test)
  kruskal_result <- kruskal.test(as.formula(paste(metabolite, "~ Group")), data = metabolite_data)
  kruskal_p_value <- kruskal_result$p.value
  kruskal_significance <- ifelse(kruskal_p_value < 0.05, "Significant", "Not Significant")
  
# Append Kruskal-Wallis test results after the redukt data frame
  result_df <- rbind(result_df, data.frame(
    Metabolite = metabolite, 
    Test = "Kruskal-Wallis", 
    p_value = kruskal_p_value, 
    Significance = kruskal_significance))
}

# Displaying the results of ANOVA and Kruskal-Wallis tests
cat("\nANOVA and Kruskal-Wallis Test Results:\n")
print(result_df)


# Print the first few lines of the data to the console
print("First few lines of the full dataset:")
print(head(metabolome_data))


# To Make sure no missing values



metabolome_data <- metabolome_data %>%
  select(Group, starts_with("Metabolite")) %>%
  na.omit()

# Make loop through each metabolite to print detailed test summaries
for (metabolite in paste0("Metabolite_", 1:10)) {
  # Extracting the metabolite data
  metabolite_data <- metabolome_data %>%
    select(Group, all_of(metabolite))
  
  # Performing ANOVA and display summary
  anova_result <- aov(as.formula(paste(metabolite, "~ Group")), data = metabolite_data)
  cat("\nANOVA summary for", metabolite, ":\n")
  print(summary(anova_result))
  
  # Performing Kruskal-Wallis test and display of results
  kruskal_result <- kruskal.test(as.formula(paste(metabolite, "~ Group")), data = metabolite_data)
  cat("\nKruskal-Wallis test result for", metabolite, ":\n")
  print(kruskal_result)
}


#Part 4: Chi-Square tests for independence

# Creating a Contingency table of Gender and Disease status
gender_disease <- table(combine_data$Gender, combine_data$Group)
chi_square_output<- chisq.test(gender_disease)

# Displaying the result
print(chi_square_output)


# Part 5: Correlation Analysis

#install.packages("ggcorrplot")
library(ggcorrplot)
#install.packages("pheatmap")
library(pheatmap) 

# Etracting species and metabolite columns for combine dataset
species_data <- combine_data %>% select(starts_with("Species"))
metabolite_data <- combine_data %>% select(starts_with("Metabolite"))

# Computimg the correlation matrix
correlation_matrix <- cor(species_data, metabolite_data, method = "pearson")
print(correlation_matrix)
correlation_plot <- ggcorrplot(correlation_matrix, lab = TRUE, lab_size = 2.5, 
                               title = "Correlation Between Microbial Species and Metabolites",
                               show.legend = TRUE)

ggsave("correlation_plot.png", plot = correlation_plot, width = 10, height = 8)

#species and metabolite columns
species_data <- combine_data %>% select(starts_with("Species_"))
metabolite_data <- combine_data %>% select(starts_with("Species_"))

# Computing the correlation matrix
#Correlation_matrix <- cor(species_data, metabolite_data, method = "pearson")

pheatmap(correlation_matrix,
         main = "Correlation Heatmap of Microbial Species and Metabolites",
         color = colorRampPalette(c("blue", "white", "yellow"))(50), # Blue to yellow gradient
         clustering_method = "average", 
         border_color = NA)


#Part 6: Multivariate Analysis (MANOVA)

#Selecting top microbial species for analysis
top_species <- combine_data %>% 
  select(Species_1, Species_2, Species_3, Species_4, Species_5)

#Selecting top metabolites for analysis
top_metabolites <- combine_data %>% 
  select(Metabolite_1, Metabolite_2, Metabolite_3, Metabolite_4, Metabolite_5)

# Combine selected species and metabolites into single dataset
independent_vars <- cbind(top_species, top_metabolites)

# Defining the dependent variables from the dataset
dependent_vars <- combine_data %>% select(BMI, Age)

# Performing the MANOVA test and the summary 
manova_result <- manova(as.matrix(dependent_vars) ~ as.matrix(independent_vars))
summary(manova_result, test = "Pillai")



